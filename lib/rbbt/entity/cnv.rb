require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/entity/gene'

Workflow.require_workflow "Sequence"

module CNV
  extend Entity
  self.annotation :jobname
  self.annotation :organism

  self.format = "Copy Number Variation"

  property :variation => :single2array do
    self.split(":").last
  end

  property :loss? => :array2single do
    @loss ||= self.variation.collect{|v| v =~/loss/i}
  end

  property :gain? => :array2single do
    @gain ||= self.variation.collect{|v| v =~/gain/i}
  end

  property :genes => :array2single do
    @genes ||= begin
                 genes = Sequence.job(:genes_at_genomic_ranges, jobname, :organism => organism, :ranges => self, :unnamed => true).run
                 genes = genes.values_at *self
                 Gene.setup(genes, "Ensembl Gene ID", organism)
               end
  end

end
 
