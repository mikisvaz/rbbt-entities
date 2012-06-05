require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/entity/gene'

Workflow.require_workflow "Sequence"

module ChromosomeRange
  extend Entity

  self.annotation :organism
 
  self.format = "Chromosome Range"


  def self.text_to_unit(text)
    text = text.sub('^', '+')
    base = text.to_f
    case
    when text =~ /KB?$/
      base * 1000
    when text =~ /MB?$/
      base * 1000_000
    when text =~ /^\d+(\.\d+)?(e\+\d+)?$/
      base
    else 
      raise "Text format not understood: #{ text }"
    end.to_i
  end

  property :unit => :array2single do
    self.collect{|range|
      chr, start, eend = range.split(":")
      [chr, text_to_unit(start), text_to_unit(eend)] * ":"
    }
  end
  persist :unit

  property :genes => :array2single do
    Sequence.job(:genes_at_genomic_ranges, "ChromosomeRange", :organism => organism, :ranges => self.unit).run.tap{|t| t.namespace = organism}.values_at *self.unit
  end

  property :chromosome => :array2single do
    self.collect{|r| r.split(":")[0] }
  end

  property :start => :array2single do
    unit.collect{|r| r.split(":")[1] }
  end

  property :eend => :array2single do
    unit.collect{|r| r.split(":")[2] }
  end

  property :ensembl_browser => :single2array do
    "http://#{Misc.ensembl_server(self.organism)}/Homo_sapiens/Location/View?db=core&r=#{chromosome}:#{start}-#{eend}"
  end
 
end
