require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')

require 'test/unit'
require 'rbbt/util/tmpfile'
require 'test/unit'
require 'rbbt/entity/gene'

class TestGene < Test::Unit::TestCase
  def test_to
    assert_equal "1020", Gene.setup("CDK5", "Associated Gene Name", "Hsa").to("Entrez Gene ID")
  end
end


