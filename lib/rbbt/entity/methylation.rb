require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/entity/chromosome_range'

module Methylation
  extend Entity
  #self.annotation :jobname
  include ChromosomeRange

  self.format = "Methylation"

  property :beta => :single2array do
    self.split(":")[3].to_f
  end

  property :methylated? => :array2single do |*args|
    threshold, *rest = args
    threshold = 0.5 if threshold.nil?
    self.beta.collect{|b| b > threshold }
  end

  property :unmethylated? => :array2single do |*args|
    threshold, *rest = args
    threshold = 0.5 if threshold.nil?
    self.beta.collect{|b| b <= threshold}
  end
end
 
