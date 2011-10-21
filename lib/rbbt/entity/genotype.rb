require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/mutation/mutation_assessor'
require 'rbbt/entity/protein'

Workflow.require_workflow "Sequence"

module MutatedIsoform
  extend Entity
  self.annotation :organism

  self.format = "Mutated Isoform"

  def protein
    if Array === self
      Protein.setup(self.collect{|mutation| mutation.split(":").first}, "Ensembl Protein ID", organism)
    else
      Protein.setup(self.split(":").first, "Ensembl Protein ID", organism)
    end
  end

  def ensembl_protein_image_url
    if Array === self
      self.collect{|e| e.ensembl_protein_image_url}
    else
      ensembl_url = if organism == "Hsa" then "www.ensembl.org" else "#{organism.sub(/.*\//,'')}.archive.ensembl.org" end
      "http://#{ensembl_url}/Homo_sapiens/Component/Transcript/Web/TranslationImage?db=core;p=#{protein};_rmd=d2a8;export=svg"
    end
  end

  ASTERISK = "*"[0]
  CONSECUENCES = %w(UTR SYNONYMOUS NOSTOP MISS-SENSE INDEL FRAMESHIFT NONSENSE)
  def single_consecuence
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

  def ary_consecuence
    self.collect{|mutation| mutation.single_consecuence}
  end


  def consecuence
    Array === self ? ary_consecuence : single_consecuence
  end

  def filter(*consecuences)
    if Array === self
      list = self.zip(consecuence).select do |mutation, consecuence|
        consecuences.include? consecuence
      end.collect{|mutation, consecuence| mutation}
    else
      list = []
      list << self if consecuences.include? consecuence
    end
    MutatedIsoform.setup(list, organism)
  end

  def self2mutation_assessor_prediction
    if Array === self
      filtered = filter "MISS-SENSE"
      correspondance = {}
      mutations = filtered.zip(filtered.protein.to "UniProt/SwissProt ID").collect do |mutation, uniprot|
        prot, change = mutation.split(":")
        next if uniprot.nil?
        uniprot_change = [uniprot, change]
        correspondance[uniprot_change] ||= []
        correspondance[uniprot_change] << mutation
        uniprot_change
      end.compact

      return TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"], :type => :list) if mutations.empty?

      tsv = MutationAssessor.chunked_predict(mutations)

      return TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"], :type => :list) if tsv.nil? or tsv.empty? 

      new = TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"], :type => :list)

      tsv.each do |key, values|
        correspondance[key.split(" ")].each do |mutation|
          new[mutation] = values
        end
     end

      new
    else
      prot, change = self.split(":")
      uniprot = protein.to "UniProt/SwissProt ID"
      return TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"], :type => :list) if uniprot.nil? or type != "MISS-SENSE"

      mutations = [[uniprot, change]]

      tsv = MutationAssessor.chunked_predict(mutations)
      tsv.add_field "Mutated Isoform" do |key, values|
        self
      end
      tsv.reorder "Mutated Isoform", ["Func. Impact"]
    end
  end

  def early_nonsense
    protein_sequences = Organism.protein_sequence(organism).tsv :persist => true, :type => :single
    filter("NONSENSE").select{|isoform_mutation|
      protein, mutation = isoform_mutation.split ":"
      if protein_sequences.include? protein
        mutation.match(/(\d+)/)[1].to_f < protein_sequences[protein].length.to_f * 0.7
      else
        Log.debug "Sequence for protein was missing: #{protein}"
        false
      end
    }
  end

  def early_frameshifts
    protein_sequences = Organism.protein_sequence(organism).tsv :persist => true, :type => :single
    filter("FRAMESHIFT").select{|isoform_mutation|
      protein, mutation = isoform_mutation.split ":"
      if protein_sequences.include? protein
        mutation.match(/(\d+)/)[1].to_f < protein_sequences[protein].length.to_f * 0.7
      else
        Log.debug "Sequence for protein was missing: #{protein}"
        false
      end
    }
  end

  def damaged(options = {})
    options = Misc.add_defaults options, :mutation_assesor_cutoff => :medium, :non_sense => true, :frameshift => true

    levels = [:low, :medium, :high].collect{|v| v.to_s}
    cutoff = levels.index options[:mutation_assesor_cutoff].to_s

    predicted = self2mutation_assessor_prediction.select{|k, v| 
      if v.nil?
        false
      else
        value = levels.index(v[0].to_s)
        value and value >= cutoff
      end
    }.collect{|k,v| k}

    predicted += early_nonsense if options[:non_sense]
    predicted += early_frameshifts if options[:frameshift]

    MutatedIsoform.setup(predicted, organism)
  end

  def damaged?(options = {})
    damaged = self.damaged
    if Array === self
      self.collect{|e| damaged.include? e}
    else
      damaged.include? self
    end
  end

  def self2damage_score
    early_nonsense = self.early_nonsense
    early_frameshifts = self.early_frameshifts
    predictions = self.self2mutation_assessor_prediction


    scores = {}
    levels = [:low, :medium, :high].collect{|v| v.to_s}
    self.each{|im|
      score = case
              when early_nonsense.include?(im)
                2
              when early_frameshifts.include?(im)
                2
              when predictions.include?(im)
                (levels.index(predictions[im].to_s) || 1)
              else 
                0
              end
      scores[im] = score
    }

    scores
  end
end

module GenomicMutation
  extend Entity
  self.annotation :jobname
  self.annotation :organism
  self.annotation :watson

  self.format = "Genomic Mutation"

  def parts
    self.split ":"
  end

  def score
    parts[3]
  end

  def self2genes
    Sequence.job(:genes_at_genomic_positions, jobname, :organism => organism, :positions => Array === self ? self : [self]).run
  end

  def genes
    Gene.setup(self2genes.values.flatten.uniq, "Ensembl Gene ID", organism)
  end

  def self2mutated_isoforms
    Sequence.job(:mutated_isoforms_for_genomic_mutations, jobname, :watson => watson, :organism => organism, :mutations => Array === self ? self : [self]).run
  end

  def self2affected_exons
    Sequence.job(:exons_at_genomic_positions, jobname, :organism => organism, :positions => Array === self ? self : [self]).run
  end

  def self2exon_junctions
    Sequence.job(:exon_junctions_at_genomic_positions, jobname, :organism => organism, :positions => Array === self ? self : [self]).run
  end

  def mutated_isoforms
    isoforms = self2mutated_isoforms.values.flatten
    MutatedIsoform.setup(isoforms, organism)
    isoforms
  end

  def in_exon_junction?
    if Array === self
      self2exon_junctions.values_at(*self).collect{|v| v.nil? ? false : v.any?}
    else
      self2exon_junctions.include?(self)? self2exon_junctions[self].any? : false
    end
  end

  def affected_exons
    self2affected_exons.values.flatten
  end

  def damaging_mutations(options = {})
    damaged_isoforms = mutated_isoforms.damaged?(options)
    damaging_mutations = self2mutated_isoforms.select{|mutation, values|
      mutated_isoforms = values
      (damaged_isoforms & mutated_isoforms).any?
    }.collect{|mutation, mutated_isoforms| mutation.dup}
    damaging_mutations + self.self2exon_junctions.reject{|mut, list| list.nil? or list.empty?}.collect{|mut, list| mut}
    GenomicMutation.setup(damaging_mutations, jobname + '.damaging', organism, watson)
  end

  def damaging?
    self.make_list.damaging_mutations.include? self
  end

  def mutations_at_genes(genes)
    genes = genes.to("Ensembl Gene ID")
    genes.compact! if Array === genes
    s2g = self.self2genes 
    subset = s2g.select("Ensembl Gene ID" => genes).keys.collect{|e| e.dup}
    GenomicMutation.setup(subset, (jobname || "Default") + '.mutations_at_genes', organism, watson)
  end

  def self2damage_score
    if Array === self
      isoform_mutation_damage = self.mutated_isoforms.self2damage_score

      damage = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Damage Score"], :type => :single, :cast => :to_f)
      self2mutated_isoforms.each{|gene, isoforms| 
        next if isoforms.empty?
        score = [isoform_mutation_damage.values_at(*isoforms).max, (gene.in_exon_junction? ? 2 : 0)].max
        damage[gene] = [damage[gene] || 0, score].max
      }
      damage
    else
      self.make_list.self2damage_score[self]
    end
  end

  def self2consecuence
    if Array === self
      self2mutated_isoforms = self.self2mutated_isoforms
      exon_junctions = self.self2exon_junctions
      consecuences = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Consecuence"], :type => :single)
      self.each do |mutation|
        mutated_isoforms = self2mutated_isoforms[mutation]
        if not mutated_isoforms.nil? and mutated_isoforms.any?
          value = mutated_isoforms.consecuence.collect{|c| MutatedIsoform::CONSECUENCES.index c}.max
          consecuences[mutation] = MutatedIsoform::CONSECUENCES[value]
        else
          consecuences[mutation] = "NONE"
        end
        consecuences[mutation] += " EXON-JUNCTION" if exon_junctions.include? mutation and exon_junctions[mutation].any?
      end
      consecuences
    else
      self.make_list.self2consecuence
    end
  end

  def consecuence
    if Array === self
      self2consecuence.values_at(*self).collect{|v| v.nil? ? "NONE" : MutatedIsoform::CONSECUENCES[v]}
    else
      self2consecuence[self] || "NONE"
    end
  end
end
