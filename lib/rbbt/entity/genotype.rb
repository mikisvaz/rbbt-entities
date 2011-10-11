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
  def single_type
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

  def ary_type
    self.collect{|mutation| mutation.single_type}
  end


  def type
    Array === self ? ary_type : single_type
  end

  def filter(*types)
    list = self.zip(type).select do |mutation, type|
      types.include? type
    end.collect{|mutation, type| mutation}

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
        correspondance[uniprot_change] = mutation
        uniprot_change
      end.compact

      tsv = MutationAssessor.chunked_predict(mutations)
      return TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"]) if tsv.nil? or tsv.empty?
      tsv.add_field "Mutated Isoform" do |key, values|
        correspondance[key.split(" ")]
      end
      tsv.reorder "Mutated Isoform", ["Func. Impact"]
    else
      prot, change = mutation.split(":")
      uniprot = protein.to "UniProt/SwissProt ID"
      mutations = [uniprot, change]

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
end

module GenomicMutation
  extend Entity
  self.annotation :name
  self.annotation :organism

  self.format = "Genomic Mutation"

  def self2genes
    Sequence.job(:genes_at_genomic_positions, name, :organism => organism, :positions => Array === self ? self : [self]).run
  end

  def genes
    Gene.setup(self2genes.values.flatten.uniq, "Ensembl Gene ID", organism)
  end

  def self2mutated_isoforms
    Sequence.job(:mutated_isoforms_for_genomic_mutations, name, :organism => organism, :mutations => Array === self ? self : [self]).run
  end

  def mutated_isoforms
    MutatedIsoform.setup(self2mutated_isoforms.values.flatten, organism)
  end

  def damaging_mutations(options = {})
    damaged_isoforms = mutated_isoforms.damaged(options)
    damaging_mutations = self2mutated_isoforms.select{|mutation, values|
      mutated_isoforms = values["Mutated Isoform"]
      (damaged_isoforms & mutated_isoforms).any?
    }.collect{|mutation, mutated_isoforms| mutation.dup}
    GenomicMutation.setup(damaging_mutations, name + '.damaging', organism)
  end

  def mutations_at_genes(genes)
    genes = genes.to("Ensembl Gene ID").compact
    s2g = self.self2genes 
    subset = s2g.select("Ensembl Gene ID" => genes).keys.collect{|e| e.dup}
    GenomicMutation.setup(subset, name + '.mutations_at_genes', organism)
  end
end
