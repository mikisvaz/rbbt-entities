require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/mutation/mutation_assessor'
require 'rbbt/entity/protein'
require 'rbbt/entity/gene'
require 'rbbt/entity/mutated_isoform'

Workflow.require_workflow "Sequence"

module GenomicMutation
  extend Entity
  self.annotation :jobname
  self.annotation :organism
  self.annotation :watson

  self.format = "Genomic Mutation"

  property :guess_watson => :array do
    if Array === self
      @watson = Sequence.job(:is_watson, jobname, :mutations => self.clean_annotations, :organism => organism).run
    else
      @watson = Sequence.job(:is_watson, jobname, :mutations => [self.clean_annotations], :organism => organism).run
    end
  end
  persist :guess_watson

  def watson
    if @watson.nil? and Array === self
      @watson = :missing
      @watson = guess_watson
    end
    @watson = false if @watson == "false"
    @watson
  end

  def orig_watson
    @watson
  end

  property :ensembl_browser => :single2array do
    "http://#{Misc.ensembl_server(self.organism)}/Homo_sapiens/Location/View?db=core&r=#{chromosome}:#{position - 100}-#{position + 100}"
  end
  persist :ensembl_browser

  property :chromosome => :array2single do
    self.clean_annotations.collect{|mut| mut.split(":")[0]}
  end
  persist :_ary_chromosome

  property :position => :array2single do
    self.clean_annotations.collect{|mut| mut.split(":")[1].to_i}
  end
  persist :_ary_position

  property :base => :array2single do
    self.clean_annotations.collect{|mut| mut.split(":")[2]}
  end
  persist :_ary_base

  property :reference => :array2single do
    Sequence.reference_allele_at_chr_positions(organism, chromosome, position)
  end
  persist :_ary_reference

  property :gene_strand_reference => :array2single do
    genes = self.genes
    gene_strand = Misc.process_to_hash(genes.compact.flatten){|list| list.strand }
    reverse = genes.collect{|list| not list.nil? and list.select{|gene| gene_strand[gene].to_s == "-1" }.any? }
    reference.zip(reverse).collect{|reference,reverse|
      reverse ? Misc::BASE2COMPLEMENT[reference] : reference
    }
  end
  persist :_ary_gene_strand_reference

  # DID NOT TRY THIS. Its supposed to deal with ambiguous gene overlaps by
  # taking the first coding gene, if any
  #property :gene_strand_reference => :array2single do
  #  genes = self.genes
  #  gene_strand = Misc.process_to_hash(genes.compact.flatten){|list| list.strand }
  #  gene_biotype = Misc.process_to_hash(genes.compact.flatten){|list| list.biotype }
  #  reference.zip(genes).collect{|reference,genes|
  #    case
  #    when (genes.nil? or genes.empty?)
  #      reference
  #    when genes.length == 1
  #      gene_strand[genes[0]] == "-1" ? Misc::BASE2COMPLEMENT[reference] : reference
  #    else
  #      coding_genes = genes.zip(gene_strand.values_at(*genes)).select{|gene,strand| gene_biotype[gene] == "protein_coding"}

  #      if coding_genes.empty?
  #        reference
  #      else
  #        coding_genes[0][1] == "-1" ? Misc::BASE2COMPLEMENT[reference] : reference
  #      end
  #    end
  #  }
  #end
  #persist :_ary_gene_strand_reference


  property :score => :array2single do
    self.clean_annotations.collect{|mut| mut.split(":")[3].to_f}
  end
  persist :_ary_score

  property :remove_score => :array2single do
    self.annotate(self.collect{|mut| mut.split(":")[0..2] * ":"})
  end
  persist :remove_score

  property :noscore => :single2array do
    self.annotate self.clean_annotations.collect{|mut| mut.split(":")[0..2]}
  end
  persist :noscore

  property :to_watson => :array2single do
    if watson
      self
    else
      result = Sequence.job(:to_watson, jobname, :mutations => self.clean_annotations, :organism => organism).run 
      self.annotate(result)
      result
    end
  end
  persist :to_watson

  property :reference => :array2single do
    Sequence.job(:reference_allele_at_genomic_positions, jobname, :positions => self.clean_annotations, :organism => organism).run.values_at *self
  end
  persist :reference

  property :type => :array2single do
    reference = watson ? self.reference : self.gene_strand_reference

    self.base.zip(reference).collect do |base,reference|

      type = case
             when base == reference
               "none"
             when (base.nil? or reference.nil? or base == "?" or reference == "?")
               "unknown"
             when (base.length > 1 or base == '-')
               "indel"
             when (not %w(A G T C).include? base and not %w(A G T C).include? reference) 
               "unknown"
             when (((Misc::IUPAC2BASE[base] || []) & ["A", "G"]).any? and     ((Misc::IUPAC2BASE[reference] || []) & ["T", "C"]).any?)
               "transversion"
             when (((Misc::IUPAC2BASE[base] || []) & ["T", "C"]).any? and     ((Misc::IUPAC2BASE[reference] || []) & ["A", "G"]).any?)
               "transversion"
             when (((Misc::IUPAC2BASE[base] || []) & ["A", "G"]).any? and     ((Misc::IUPAC2BASE[reference] || [nil]) & ["T", "C", nil]).empty?)
               "transition"
             when (((Misc::IUPAC2BASE[base] || []) & ["T", "C"]).any? and     ((Misc::IUPAC2BASE[reference] || [nil]) & ["A", "G", nil]).empty?)
               "transition"
             else
               "unknown"
             end
      type
    end

  end
  persist :type

  property :offset_in_genes => :array2single do
    gene2chr_start = Misc.process_to_hash(genes.flatten){|list| list.chr_start}
    position.zip(genes).collect{|position, list|
      list.collect{|gene|
        next if not gene2chr_start.include? gene
        [gene, position.to_i - gene2chr_start[gene]] * ":"
      }.compact
    }
  end
  persist :offset_in_genes

  property :genes => :array2single do
    genes = Sequence.job(:genes_at_genomic_positions, jobname, :organism => organism, :positions => self.clean_annotations).run
    genes.unnamed = true
    genes = genes.values_at *self
    Gene.setup(genes, "Ensembl Gene ID", organism)
  end
  persist :genes

  property :mutated_isoforms => :array2single do
    res = Sequence.job(:mutated_isoforms_for_genomic_mutations, jobname, :watson => watson, :organism => organism, :mutations => self.clean_annotations).run.values_at *self
    res.each{|list| list.organism = organism unless list.nil?}
    res.compact[0].annotate res if res.compact[0].respond_to? :annotate
    res
  end
  persist :mutated_isoforms

  property :exon_junctions => :array do
    Sequence.job(:exon_junctions_at_genomic_positions, jobname, :organism => organism, :positions => self.clean_annotations).run.values_at(*self)
  end
  persist :exon_junctions

  property :in_exon_junction? => :array2single do
    exon_junctions.collect{|l| not l.nil? and not l.empty?}
  end
  persist :in_exon_junction?

  property :over_range? => :array2single do |range|
    chromosome.zip(position).collect{|chr,pos| chr == gene_chromosome and range.include? pos}
  end
  persist :_ary_over_range?

  property :over_gene? => :array2single do |gene|
    gene = Gene.setup(gene.dup, "Ensembl Gene ID", organism) unless Gene === gene

    gene_range = gene.chr_range
    gene_chromosome = gene.chromosome

    if gene_range.nil?
      [false] * self.length
    else
      chromosome.zip(position).collect{|chr,pos| chr == gene_chromosome and gene_range.include? pos}
    end
  end
  persist :_ary_over_gene?

  property :affected_exons  => :array2single do
    Sequence.job(:exons_at_genomic_positions, jobname, :organism => organism, :positions => self.clean_annotations).run.values_at *self
  end
  persist :affected_exons

  property :damaging? => :array2single do |*args|
    damaged_mutated_isoforms = mutated_isoforms.compact.flatten.select{|mi| mi.damaged?(*args)}
    exon_junctions.zip(mutated_isoforms).collect do |exs, mis|
      (Array === exs and exs.any?) or
      (Array === mis and (damaged_mutated_isoforms & mis).any?)
    end
  end
  persist :damaging?
end
