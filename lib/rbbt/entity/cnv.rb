require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/entity/chromosome_range'

module CNV
  extend Entity
  #self.annotation :jobname
  include ChromosomeRange

  self.format = "Copy Number Variation"

  property :variation => :single2array do
    var = self.split(":")[3]
    if var =~ /^-?\d+$/
      case
      when var.to_i <= -1 
        "loss"
      when var.to_i >= 1 
        "gain"
      else
        ""
      end
    else
      var
    end
  end

  property :loss? => :array2single do
    self.variation.collect{|v| v =~/loss/i}
  end

  property :gain? => :array2single do
    self.variation.collect{|v| v =~/gain/i}
  end

  #property :genes => :array2single do
  #  @genes ||= begin
  #               genes = Sequence.job(:genes_at_genomic_ranges, jobname, :organism => organism, :ranges => self, :unnamed => true).run
  #               genes = genes.chunked_values_at self
  #               Gene.setup(genes, "Ensembl Gene ID", organism)
  #             end
  #end

  #property :chromosome => :array2single do
  #  self.clean_annotations.collect{|mut| mut.split(":")[0]}
  #end

  #property :start => :array2single do
  #  self.clean_annotations.collect{|mut| mut.split(":")[1].to_i}
  #end

  #property :begin => :array2single do
  #  self.start
  #end

  #property :end => :array2single do
  #  self.clean_annotations.collect{|mut| mut.split(":")[2].to_i}
  #end

  #property :range => :array2single do
  #  start.zip(self.end).collect{|s,e| (s..e)}
  #end
end
 
