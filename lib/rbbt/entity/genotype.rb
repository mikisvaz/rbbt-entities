require 'rbbt/workflow'
require 'rbbt/entity'
require 'rbbt/entity/genomic_mutation'

module Genotype
  extend Workflow

  if self.respond_to? :extended
    class << self
      alias prev_genotype_extended extended
    end
  end

  def self.extended(base)
    prev_genotype_extended(base) if self.respond_to? :prev_genotype_extended
    base.helper :genotype do
      base
    end
  end

  module Cohort
    extend Workflow

    if self.respond_to? :extended
      class << self
        alias prev_genotype_cohort_extended extended
      end
    end

    def self.extended(cohort)
      prev_genotype_cohort_extended(cohort) if self.respond_to? :prev_genotype_cohort_extended

      class << cohort
        attr_accessor :metagenotype

        def jobname
          if @jobname.nil?
            @jobname ||= "Meta-genotype: " + self.collect{|g| g.jobname} * ", "
            @jobname[100..-1] = " (etc; #{self.length} genotypes)" if @jobname.length > 100
          end
          @jobname
        end

        def metagenotype
          if @metagenotype.nil?
            @metagenotype = GenomicMutation.setup(self.dup.flatten, jobname, self[0].organism, self[0].orig_watson)
            @metagenotype.extend Genotype unless Genotype === @metagenotype
          end
          @metagenotype
        end
      end unless cohort.respond_to? :metagenotype

      cohort.each do |genotype| genotype.extend Genotype unless Genotype === genotype end

      cohort.helper :metagenotype do
        cohort.metagenotype
      end

      cohort.helper :samples do
        cohort
      end

      NamedArray.setup(cohort, cohort.collect{|genotype| genotype.jobname})
    end

    def subset(genotypes)
      new = self.values_at *(genotypes & fields)
      new.extend Cohort
    end

    returns "Ensembl Gene ID"
    task :all_affected_genes => :array do
      set_info :organism, metagenotype.organism
      samples.collect{|genotype| genotype.all_affected_genes}.flatten.uniq
    end

    returns "Ensembl Gene ID"
    input :methods, :array, "Predictive methods", [:sift, :mutation_assessor]
    input :threshold, :float, "from 0 to 1", 0.8
    task :damaged_genes => :array do |methods, threshold|
      set_info :organism, metagenotype.organism
      samples.collect{|genotype| genotype.damaged_genes(:methods => methods, :threshold => threshold)}.flatten.uniq
    end

    returns "Ensembl Gene ID"
    task :recurrent_genes => :array do
      set_info :organism, metagenotype.organism
      count = Hash.new(0)
      samples.each do |genotype| genotype.genes.flatten.uniq.each{|gene| count[gene] += 1} end
      count.select{|gene, c| c > 1}.collect{|gene,c| gene.dup}
    end

    %w(damaged_genes recurrent_genes all_affected_genes).each do |name|
      define_method name do |*args|
        options = args.first
        @cache ||= {}
        key = [name, Misc.hash2md5(options || {})]
        @cache[key] ||= self.job(name, self.jobname, options || {}).run
      end
    end

  end

  returns "Ensembl Gene ID"
  task :all_affected_genes => :array do 
    set_info :organism, genotype.organism
    genotype.genes.clean_annotations.flatten.uniq
  end

  dep :all_affected_genes
  returns "Ensembl Gene ID"
  task :long_genes => :array do
    all_affected_genes = step(:all_affected_genes).load
    long_genes = all_affected_genes.select{|gene| 
      length = gene.max_protein_length
      length and length > 1000 or gene.name =~ /^PCDH/
    }

    set_info :organism, genotype.organism
    long_genes.clean_annotations
  end

  returns "Ensembl Gene ID"
  task :mutations_in_exon_junctions => :array do
    set_info :organism, genotype.organism
    genotype.select{|mutation| mutation.in_exon_junction?}.clean_annotations
  end

  returns "Ensembl Gene ID"
  task :with_non_synonymous_mutations => :array do
    set_info :organism, genotype.organism
    genotype.mutated_isoforms.flatten.compact.reject{|mutated_isoform| ["SYNONYMOUS", "UTR"].include? mutated_isoform.consequence}.transcript.gene.uniq
  end

  returns "Ensembl Gene ID"
  input :methods, :array, "Predictive methods", [:sift, :mutation_assessor]
  input :threshold, :float, "from 0 to 1", 0.8
  task :with_damaged_isoforms => :array do |methods,threshold|
    set_info :organism, genotype.organism
    mutated_isoform_damage = Misc.process_to_hash(genotype.mutated_isoforms.flatten.compact){|list| MutatedIsoform.setup(list, genotype.organism).damage_scores(methods)}
    genotype.select{|mutation|  if mutation.mutated_isoforms then mutated_isoform_damage.values_at(*mutation.mutated_isoforms.flatten.compact).select{|score| not score.nil?  and score > threshold}.any? else false; end}.genes.flatten.uniq.clean_annotations
  end

  returns "Ensembl Gene ID"
  task :truncated => :array do
    set_info :organism, genotype.organism
    truncated_isoforms = MutatedIsoform.setup(genotype.mutated_isoforms.flatten.compact, "Hsa/jun2011").select{|isoform_mutation| isoform_mutation.truncated }
    proteins = truncated_isoforms.protein
    genes = proteins.gene
    genes.to("Ensembl Gene ID").uniq.clean_annotations
  end

  returns "Ensembl Gene ID"
  task :affected_exon_junctions => :array do
    set_info :organism, genotype.organism
    genotype.select{|mutation| mutation.in_exon_junction?}.genes.flatten.clean_annotations
  end

  dep :with_damaged_isoforms, :truncated, :affected_exon_junctions
  returns "Ensembl Gene ID"
  task :damaged_genes => :array do
    set_info :organism, genotype.organism

    with_damaged_isoforms   = step(:with_damaged_isoforms).load.clean_annotations
    truncated               = step(:truncated).load.clean_annotations
    affected_exon_junctions = step(:affected_exon_junctions).load.clean_annotations

    (with_damaged_isoforms + truncated + affected_exon_junctions).uniq
  end

  %w(all_affected_genes damaged_genes truncated with_damaged_isoforms with_non_synonymous_mutations affected_exon_junctions long_genes recurrent_genes).each do |name|
    define_method name do |*args|
      options = args.first
      @cache ||= {}
      key = [name, Misc.hash2md5(options || {})]
      @cache[key] ||= self.job(name, self.jobname, options || {}).run
    end
  end
end

