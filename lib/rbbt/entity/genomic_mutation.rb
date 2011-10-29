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

  property :genes => :array2single do
    @genes ||= begin
                 genes = Sequence.job(:genes_at_genomic_positions, jobname, :organism => organism, :positions => self).run
                 genes.unnamed = true
                 genes = genes.values_at *self
                 Gene.setup(genes, "Ensembl Gene ID", organism)
               end
  end

  property :mutated_isoforms => :array2single do
    @mutated_isoforms ||= begin
                            res = Sequence.job(:mutated_isoforms_for_genomic_mutations, jobname, :watson => watson, :organism => organism, :mutations => self).run.values_at *self
                            res.each{|list| list.organism = organism}
                            res
                          end
  end

  property :exon_junctions do
    @exon_junctions ||= Sequence.job(:exon_junctions_at_genomic_positions, jobname, :organism => organism, :positions => self).run.values_at *self
  end

  property :in_exon_junction? => :array2single do
    exon_junctions.collect{|l| not l.nil? and not l.empty?}
  end

  property :over_gene? => :array2single do |gene|
    @over_genes ||= {}
    @over_genes[gene] ||= genes.clean_annotations.collect{|list| list.include? gene}
  end

  property :mutation_assessor_scores => :array2single do
    @mutation_assessor_scores ||= begin
                                   mutated_isoforms = self.mutated_isoforms
                                   all_mutated_isoforms = MutatedIsoform.setup(mutated_isoforms.flatten.compact, organism)
                                   mutated_isoform2damage_score = Misc.process_to_hash(all_mutated_isoforms){|list| all_mutated_isoforms.mutation_assessor_scores}

                                   mutated_isoforms.collect{|list| list.nil? ? [] : mutated_isoform2damage_score.values_at(*list)}
                                 end
  end

  property :truncated do
    @truncated ||= begin
                     mutated_isoforms = self.mutated_isoforms
                     all_mutated_isoforms = MutatedIsoform.setup(mutated_isoforms.flatten.compact, organism)
                     mutated_isoform2truncated = Misc.process_to_hash(all_mutated_isoforms){|list| all_mutated_isoforms.truncated}
                     mutated_isoforms.collect{|list| list.nil? ? [] : mutated_isoform2truncated.values_at(*list)}
                   end
  end


end
