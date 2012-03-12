require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/mutation/mutation_assessor'
require 'rbbt/mutation/sift'
require 'rbbt/entity/protein'
require 'rbbt/sources/uniprot'
require 'rbbt/entity/gene'
require 'nokogiri'

Workflow.require_workflow 'structure'

module MutatedIsoform
  extend Entity
  self.annotation :organism

  self.format = "Mutated Isoform"

  property :protein => :array2single do
    Protein.setup(self.collect{|mutation| mutation.split(":").first if mutation =~ /^ENSP/}, "Ensembl Protein ID", organism)
  end
  persist :protein

  property :transcript => :array2single do
    begin
      protein = self.protein
      Transcript.setup(protein.transcript.zip(self.collect{|mutation| mutation.split(":").first}).collect{|p| p.compact.first}, "Ensembl Transcript ID", organism)
    end
  end
  persist :transcript

  property :change => :array2single do
    self.collect{|mi| mi.split(":").last}
  end
  persist :change

  property :position => :array2single do
    change.collect{|c|
      if c.match(/[^\d](\d+)[^\d]/)
        $1.to_i
      else
        nil
      end
    }
  end
  persist :position

  property :ensembl_protein_image_url => :single2array do
    ensembl_url = if organism == "Hsa" then "www.ensembl.org" else "#{organism.sub(/.*\//,'')}.archive.ensembl.org" end
    "http://#{ensembl_url}/Homo_sapiens/Component/Transcript/Web/TranslationImage?db=core;p=#{protein};_rmd=d2a8;export=svg"
  end
  persist :ensembl_protein_image_url

  property :marked_svg => :single2array do
    svg = Open.read(protein.ensembl_protein_image_url)
    
    seq_len = protein.sequence_length
    position = self.position

    doc = Nokogiri::XML(svg)
    return nil unless doc.css('svg')
    width = doc.css('svg').first.attr('width').to_f
    height = doc.css('svg').first.attr('height').to_f
    start = doc.css('rect.ac').first.attr('x').to_f

    if width and height and start and seq_len and position
      offset = (width - start)/seq_len * position + start
      svg.sub(/<\/svg>/,"<rect x='#{offset}' y='1' width='1' height='#{height}' style='fill:rgb(255,0,0);opacity:0.5;stroke:none;'></svg>")
    else
      svg
    end
  end
  persist :marked_svg

  ASTERISK = "*"[0]
  CONSECUENCES = %w(UTR SYNONYMOUS NOSTOP MISS-SENSE INDEL FRAMESHIFT NONSENSE)
  property :consequence => :single2array do
    prot, change = self.split(":")

    case
    when change =~ /UTR/
      "UTR"
    when (change[0] == ASTERISK and not change[0] == change[-1])
      "NOSTOP"
    when (change[-1] == ASTERISK and not change[0] == change[-1])
      "NONSENSE"
    when change =~ /Indel/
      "INDEL"
    when change =~ /FrameShift/
      "FRAMESHIFT"
    when change[0] == change[-1]
      "SYNONYMOUS"
    else
      "MISS-SENSE"
    end
  end
  persist :consequence

  property :truncated => :array2single do
    begin
      proteins = self.protein.compact.flatten
      protein2sequence_length = Misc.process_to_hash(proteins){|list| proteins.any? ? proteins.sequence_length : []}

      self.collect do |isoform_mutation|

                       next if isoform_mutation.consequence != "FRAMESHIFT" and isoform_mutation.consequence != "NONSENSE"
                       protein  = isoform_mutation.protein
                       position = isoform_mutation.position
                       sequence_length = protein2sequence_length[protein]

                       case
                       when (sequence_length.nil? or position.nil?)
                         nil
                       when position < sequence_length.to_f * 0.7
                         true
                       else
                         false
                       end
                     end
                   end
  end
  persist :truncated

  property :damage_scores => :array2single do |*args|
    begin
      methods = args.first
      methods = [:sift, :mutation_assessor] if methods.nil?
      methods = [methods] unless Array === methods
      values = methods.collect{|method|
        case method.to_sym
        when :sift
          sift_scores
        when :mutation_assessor
          mutation_assessor_scores
        else
          raise "Unknown predictive method: #{ method }"
        end
      }
      if values.compact.empty?
        return [nil] * self.length
      else
        scores = values.shift
        scores = scores.zip(*values)

        scores.collect{|p|
          p = p.compact
          if p.empty?
            nil
          else
            p.inject(0.0){|acc, e| acc += e} / p.length
          end
        }
      end
    end
  end
  persist :damage_scores

  property :damaged? => :array2single do |*args|
    begin
      methods, threshold = args
      threshold     = 0.8 if threshold.nil?
      damage_scores = self.damage_scores(methods)
      truncated     = self.truncated
      damage_scores.zip(truncated).collect{|damage, truncated| truncated or (not damage.nil? and damage > threshold) }
    end
  end
  persist :damaged?

  property :sift_scores => :array2single do
    begin
      missense = self.select{|iso_mut| iso_mut.consequence == "MISS-SENSE"}

      values = SIFT.chunked_predict(missense).values_at(*self).collect{|v|
        v.nil? ? nil : 1.0 - v["Score 1"].to_f
      }

      values

      #range = {nil => nil,
      #  ""  => nil,
      #  "TOLERATED" => 0,
      #  "*DAMAGING" => 1,
      #  "DAMAGING" => 1}

      #range.values_at *values
    end
  end
  persist :sift_scores

  property :mutation_assessor_scores => :array2single do
    begin
      missense = self.select{|mutation| mutation.consequence == "MISS-SENSE"}

      correspondance = {}
      list = missense.zip(missense.protein.to "UniProt/SwissProt ID").collect do |mutation, uniprot|
        prot, change = mutation.split(":")
        next if uniprot.nil?
        uniprot_change = [uniprot.upcase, change.upcase]
        correspondance[uniprot_change] ||= []
        correspondance[uniprot_change] << mutation
        uniprot_change
      end.compact

      return [nil] * self.length if list.empty?

      tsv = MutationAssessor.chunked_predict(list.sort_by{|p| p * "_"})

      return [nil] * self.length if tsv.empty?

      new = TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"], :type => :list)

      tsv.each do |key, values|
        uniprot, change = key.split(" ")
        uniprot_change = [uniprot.upcase, change.upcase]
        
        if correspondance.include? uniprot_change
          correspondance[uniprot_change].each do |mutation|
            new[mutation] = values["Func. Impact"]
          end
        else
          Log.medium "Correspondace value missing: #{uniprot_change.inspect}"
        end
      end


      range = {nil => nil,
        ""  => nil,
        "neutral" => 0,
        "low" => 0.5,
        "medium" => 0.7,
        "high" => 1.0}

      range.values_at *new.values_at(*self)
    end
  end
  persist :mutation_assessor_scores

  property :pdbs => :single do
    uniprot = self.transcript.protein.uniprot
    next if uniprot.nil?
    UniProt.pdbs_covering_aa_position(uniprot, self.position)
  end
  persist :pdbs

  property :pdbs_and_positions => :single do
    pdbs.collect do |pdb, info|
      [pdb, Structure.job(:sequence_position_in_pdb, "Protein: #{ self }", :sequence => protein.sequence, :organism => organism, :position => position, :pdb => pdb).run]
    end
  end
end

