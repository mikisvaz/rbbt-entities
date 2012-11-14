require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/mutation/mutation_assessor'
require 'rbbt/mutation/sift'
require 'rbbt/entity/protein'
require 'rbbt/sources/uniprot'
require 'rbbt/sources/InterPro'
require 'rbbt/entity/gene'
require 'nokogiri'

Workflow.require_workflow 'structure'
Workflow.require_workflow 'MutEval'

module MutatedIsoform
  extend Entity
  self.annotation :organism

  self.format = "Mutated Isoform"

  DEFAULT_DAMAGE_PREDICTORS = [:sift, :mutation_assessor]

  property :protein => :array2single do
    proteins = self.collect{|mutation| mutation.split(":").first if mutation[0..3] == "ENSP"}
    Protein.setup(proteins, "Ensembl Protein ID", organism)
  end

  property :transcript => :array2single do
    begin
      protein = self.protein
      Transcript.setup(protein.transcript.zip(self.collect{|mutation| mutation.split(":").first}).collect{|p| p.compact.first}, "Ensembl Transcript ID", organism)
    end
  end

  property :change => :array2single do
    self.collect{|mi| mi.split(":").last}
  end

  property :position => :array2single do
    change.collect{|c|
      if c.match(/[^\d](\d+)[^\d]/)
        $1.to_i
      else
        nil
      end
    }
  end

  property :ensembl_protein_image_url => :single2array do
    ensembl_url = if organism == "Hsa" then "www.ensembl.org" else "#{organism.sub(/.*\//,'')}.archive.ensembl.org" end
    "http://#{ensembl_url}/Homo_sapiens/Component/Transcript/Web/TranslationImage?db=core;p=#{protein};_rmd=d2a8;export=svg"
  end

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

  ASTERISK = "*"[0]
  CONSECUENCES = %w(UTR SYNONYMOUS NOSTOP MISS-SENSE INDEL FRAMESHIFT NONSENSE)
  property :consequence => :single2array do
    return nil if self.nil?

    prot, change = self.split(":")

    case
    when change.nil?
      nil
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

  property :in_utr => :array2single do
    consequence.collect{|c|
      c == "UTR"
    }
  end

  property :synonymous => :array2single do
    consequence.collect{|c|
      c == "SYNONYMOUS"
    }
  end

  property :non_synonymous => :array2single do
    consequence.collect{|c|
      not c.nil? and c != "SYNONYMOUS" and c != "UTR"
    }
  end

  property :affected_interpro_domains => :single do
    if protein.nil?
      []
    else
      InterProDomain.setup(Misc.zip_fields(protein.interpro_domain_positions || []).select{|d,s,e|
        e.to_i > position and s.to_i < position
      }.collect{|d,s,e| d }, organism)
    end
  end

  property :affected_interpro_domain_positions => :single do
    if protein.nil?
      []
    else
      Misc.zip_fields(protein.interpro_domain_positions || []).select{|d,s,e|
        e.to_i > position and s.to_i < position
      }.collect{|d,s,e| [d, position - s.to_i, s.to_i, e.to_i]}
    end
  end

  property :affected_domain_positions => :single do
    affected_interpro_domain_positions
  end

  property :affected_domains => :single do
    affected_interpro_domains
  end

  property :ablated_interpro_domains => :single do
    if protein.nil?
      []
    else
      InterProDomain.setup(Misc.zip_fields(protein.interpro_domain_positions || []).select{|d,s,e|
        e.to_i > position
      }.collect{|d,s,e| d }, organism)
    end
  end

  property :ablated_interpro_domain_positions => :single do
    if protein.nil?
      []
    else
      Misc.zip_fields(protein.interpro_domain_positions || []).select{|d,s,e|
        e.to_i > position
      }.collect{|d,s,e| [d, s.to_i, e.to_i]}
    end
  end

  property :ablated_domain_positions => :single do
    ablated_interpro_domain_positions
  end

  property :ablated_domains => :single do
    ablated_interpro_domains
  end

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
        when (isoform_mutation.ablated_domains.any?)
          true
        else
          false
        end
      end
    end
  end

  property :damage_scores => :array2single do |*args|
    begin
      methods = args.first
      #methods = [:sift, :mutation_assessor] if methods.nil?
      methods = MutatedIsoform::DEFAULT_DAMAGE_PREDICTORS if methods.nil?
      methods = [methods] unless Array === methods
      values = methods.collect{|method|
        case method.to_sym
        when :sift
          sift_scores
        when :mutation_assessor
          mutation_assessor_scores
        when :polyphen
          polyphen_scores
        when :snps_and_go
          snps_and_go_scores
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

  property :damaged? => :array2single do |*args|
    begin
      methods, threshold = args
      threshold, methods = methods, nil if threshold.nil? and not Array === methods
      threshold     = 0.8 if threshold.nil?
      damage_scores = self.damage_scores(methods)
      truncated     = self.truncated

      damage_scores.zip(truncated).collect{|damage, truncated| truncated or (not damage.nil? and damage > threshold) }
    end
  end

  property :snps_and_go_scores => :array2single do
    begin
      missense = self.select{|mutation| mutation.consequence == "MISS-SENSE"}
      res = MutEval.job(:snps_and_go, "MutatedIsoform", :mutations => missense.sort, :organism => organism).run
      res.values_at(*self).collect{|v| (v.nil? or v["SNPSandGO Score"].nil? or v["SNPSandGO Score"].empty?) ? 
        nil : 
        (v["SNPSandGO Prediction"] == "Disease" ? 1.0 - (10.0 - v["SNPSandGO Score"].to_f) / 20 : 0 + (10.0 - v["SNPSandGO Score"].to_f) / 20)
      }
    rescue
      Log.warn $!.message
      [nil] * self.length
    end
  end

  property :polyphen_scores => :array2single do
    begin
      missense = self.select{|mutation| mutation.consequence == "MISS-SENSE"}
      res = MutEval.job(:polyphen, "MutatedIsoform", :mutations => missense.sort, :organism => organism).run
      res.values_at(*self).collect{|v| (v.nil? or v["Polyphen Score"].nil? or v["Polyphen Score"].empty?) ? nil : v["Polyphen Score"].to_f / 10}
    rescue
      Log.warn $!.message
      [nil] * self.length
    end
  end

  property :sift_scores => :array2single do
    begin
      missense = self.select{|mutation| mutation.consequence == "MISS-SENSE"}
      res = MutEval.job(:sift, "MutatedIsoform", :mutations => missense.sort, :organism => organism).run
      res.values_at(*self).collect{|v| (v.nil? or v["SIFT Score"].nil? or v["SIFT Score"].empty?) ? nil : 1.0 - v["SIFT Score"].to_f}
    rescue
      Log.warn $!.message
      [nil] * self.length
    end
  end

  property :mutation_assessor_scores => :array2single do
    range = {nil => nil,
      ""  => nil,
      "neutral" => 0,
      "low" => 0.5,
      "medium" => 0.7,
      "high" => 1.0}

    begin
      missense = self.select{|mutation| mutation.consequence == "MISS-SENSE"}
      MutEval.job(:mutation_assessor, "MutatedIsoform", :mutations => missense.sort, :organism => organism).run.values_at(*self).collect{|v| (v.nil? or v["Mutation Assessor Prediction"].nil? or v["Mutation Assessor Prediction"].empty?) ? nil : range[v["Mutation Assessor Prediction"]]}
    rescue
      Log.warn $!.message
      [nil] * self.length
    end
  end

  property :pdbs => :single do
    uniprot = self.transcript.protein.uniprot
    next if uniprot.nil?
    UniProt.pdbs_covering_aa_position(uniprot, self.position)
  end

  property :pdbs_and_positions => :single do
    pdbs.collect do |pdb, info|
      [pdb, Structure.job(:sequence_position_in_pdb, "Protein: #{ self }", :sequence => protein.sequence, :organism => organism, :position => position, :pdb => pdb).run]
    end
  end
end

