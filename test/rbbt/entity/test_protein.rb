require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')

require 'test/unit'
require 'rbbt/util/tmpfile'
require 'test/unit'
require 'rbbt/entity/protein'
require 'rbbt/entity/mutated_isoform'

class TestProtein < Test::Unit::TestCase
  PROTEIN = Protein.setup("ENSP00000275493", "Hsa/jun2011")
  PROTEIN_ARRAY = Protein.setup(["ENSP00000275493"], "Hsa/jun2011")

  def test_clean_annotations
    assert Protein === PROTEIN
    assert(!(Protein === PROTEIN.clean_annotations))
    assert Gene === PROTEIN.gene
    assert(!(Protein === PROTEIN.gene))
    assert(PROTEIN_ARRAY.respond_to? :annotated_array_clean_each)
    assert(!(PROTEIN_ARRAY.clean_annotations.respond_to? :annotated_array_clean_each))
  end

  def test_gene
    assert_equal "EGFR", PROTEIN.gene.to("Associated Gene Name")
  end

  def test_protein
    ddd MutatedIsoform.setup(["ENSP00000322422:K168FrameShift"], "Hsa/may2009").info
    ddd MutatedIsoform.setup(["ENSP00000322422:K168FrameShift"], "Hsa/may2009").protein
  end

end


