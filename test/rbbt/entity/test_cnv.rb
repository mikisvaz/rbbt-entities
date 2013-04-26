require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')
require 'test/unit'
require 'rbbt'
require 'rbbt/entity'
require 'rbbt/entity/cnv'


class TestCNV < Test::Unit::TestCase

  def test_range
    assert_equal (100..200), CNV.setup("1:100:200:-1", "Hsa/jun2011").range
  end

  def test_genes
    assert  CNV.setup("1:100MB:200MB:-1", "Hsa/jun2011").unit.genes.name.length > 0
    #assert_equal [], CNV.setup("1:100MB:200MB:-1", "Hsa/jun2011").unit.genes
  end
end

