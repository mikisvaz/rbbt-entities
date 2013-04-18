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
                 genes = genes.chunked_values_at self
                 Gene.setup(genes, "Ensembl Gene ID", organism)
               end
  end

  property :chromosome => :array2single do
    self.clean_annotations.collect{|mut| mut.split(":")[0]}
  end

  property :start => :array2single do
    self.clean_annotations.collect{|mut| mut.split(":")[1].to_i}
  end

  property :begin => :array2single do
    self.start
  end


  property :end => :array2single do
    self.clean_annotations.collect{|mut| mut.split(":")[2].to_i}
  end

  property :range => :array2single do
    start.zip(self.end).collect{|s,e| (s..e)}
  end

  property :over_chromosome_range? => :array2single do |r|
    range_chr, start, eend = r.split(":")
    range_range = (start.to_i..eend.to_i)
    chromosome.zip(range).collect do |chr,range| 
      chr == range_chr and ( range.include?(range_range.begin) or range.include?(range_range.end) or range_range.include?(range.begin) or range_range.include?(range.end))
    end
  end
end
 
