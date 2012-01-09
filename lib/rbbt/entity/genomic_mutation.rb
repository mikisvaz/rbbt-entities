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

  def watson
    if @watson.nil?
      if Array === self
        @watson = Sequence.job(:is_watson, jobname, :mutations => self.clean_annotations, :organism => organism).run
      else
        @watson = Sequence.job(:is_watson, jobname, :mutations => [self.clean_annotations], :organism => organism).run
      end
    end
    @watson
  end

  property :ensembl_browser => :single2array do
    "http://#{Misc.ensembl_server(self.organism)}/Homo_sapiens/Location/View?db=core&r=#{chromosome}:#{position - 100}-#{position + 100}"
  end

  property :chromosome => :array2single do
    @chromosome ||= self.clean_annotations.collect{|mut| mut.split(":")[0]}
  end

  property :position => :array2single do
   @position ||= self.clean_annotations.collect{|mut| mut.split(":")[1].to_i}
  end

  property :base => :array2single do
    @base ||= self.clean_annotations.collect{|mut| mut.split(":")[2]}
  end

  property :reference => :array2single do
    @reference ||= Sequence.reference_allele_at_chr_positions(organism, chromosome, position)
  end

  property :score => :array2single do
    @base ||= self.clean_annotations.collect{|mut| mut.split(":")[3].to_f}
  end

  property :remove_score => :array2single do
    @remove_score ||= self.annotate(self.collect{|mut| mut.split(":")[0..2] * ":"})
  end

  property :noscore => :single2array do
    self.annotate self.clean_annotations.collect{|mut| mut.split(":")[0..2]}
  end

  property :to_watson => :array2single do
    if watson
      self
    else
      result = Sequence.job(:to_watson, jobname, :mutations => self.clean_annotations, :organism => organism).run 
      self.annotate(result)
      result
    end
  end

  property :reference => :array2single do
    @reference ||= Sequence.job(:reference_allele_at_genomic_positions, jobname, :positions => self.clean_annotations, :organism => organism).run.values_at *self
  end

  property :type => :array2single do
    self.to_watson.base.zip(reference).collect do |base,reference|
      case
      when base == reference
        "none"
      when (base.length > 1 or base == '-')
        "indel"
      when (not %w(A G T C).include? base and not %w(A G T C).include? reference) 
        nil
      when ((base == "A" or base == "G") and (reference == "A" or reference == "G"))
        "transition"
      when ((base == "T" or base == "C") and (reference == "T" or reference == "C"))
        "transition"
      else
        "transversion"
      end
    end
  end

  property :offset_in_genes => :array2single do
    gene2chr_start = Misc.process_to_hash(genes.flatten){|list| list.chr_start}
    position.zip(genes).collect{|position, list|
      list.collect{|gene|
        next if not gene2chr_start.include? gene
        [gene, position.to_i - gene2chr_start[gene]] * ":"
      }.compact
    }
  end

  property :genes => :array2single do
    @genes ||= begin
                 genes = Sequence.job(:genes_at_genomic_positions, jobname, :organism => organism, :positions => self.clean_annotations).run
                 genes.unnamed = true
                 genes = genes.values_at *self
                 Gene.setup(genes, "Ensembl Gene ID", organism)
               end
  end

  property :mutated_isoforms => :array2single do
    @mutated_isoforms ||= begin
                            res = Sequence.job(:mutated_isoforms_for_genomic_mutations, jobname, :watson => watson, :organism => organism, :mutations => self.clean_annotations).run.values_at *self
                            res.each{|list| list.organism = organism unless list.nil?}
                            res[0].annotate res if res[0].respond_to? :annotate
                            res
                          end
  end

  property :exon_junctions do
    @exon_junctions ||= Sequence.job(:exon_junctions_at_genomic_positions, jobname, :organism => organism, :positions => self.clean_annotations).run.values_at *self
  end

  property :in_exon_junction? => :array2single do
    exon_junctions.collect{|l| not l.nil? and not l.empty?}
  end

  property :over_gene? => :array2single do |gene|
    @over_genes ||= {}
    @over_genes[gene] ||= begin
                            gene = gene.ensembl if gene.respond_to? :ensembl
                            genes.clean_annotations.collect{|list| list.include? gene}
                          end
  end

  property :affected_exons  => :array2single do
    @affected_exons ||= begin
                          Sequence.job(:exons_at_genomic_positions, jobname, :organism => organism, :positions => self.clean_annotations).run.values_at *self
                        end
  end

  property :damaging? => :array2single do
    @damaging ||= begin
                    exon_junctions.zip(mutated_isoforms).collect do |exs, mis|
                      (Array === exs and exs.any?) or
                      (Array === mis and mis.select{|mi| mi.damaged?([:sift])}.any?)
                    end
                  end
  end
end
