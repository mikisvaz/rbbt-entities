require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')

require 'test/unit'
require 'rbbt/util/tmpfile'
require 'test/unit'
require 'rbbt/entity/mutated_isoform'

class TestMutatedIsoform < Test::Unit::TestCase
  MUTATION = MutatedIsoform.setup("ENSP00000275493:G719A", "Hsa/jun2011")

  def _test_protein
    assert_equal "EGFR", MUTATION.protein.gene.to("Associated Gene Name")
  end

  def test_truncated

    protein = Protein.setup(Gene.setup("CDK5", "Associated Gene Name", "Hsa/jun2011").to("Ensembl Protein ID"), "Ensembl Protein ID", "Hsa/jun2011")

    change_position = (protein.sequence_length.to_f * 0.5).to_i
    wildtype = protein.sequence[(change_position..change_position)]
    mutation = "*"
    new_mutation = wildtype << change_position.to_s << mutation
    mutation = MutatedIsoform.setup([protein, new_mutation] * ":", "Hsa/jun2011")
    assert mutation.truncated


    change_position = (protein.sequence_length.to_f * 0.9).to_i
    wildtype = protein.sequence[(change_position..change_position)]
    mutation = "*"
    new_mutation = wildtype << change_position.to_s << mutation
    mutation = MutatedIsoform.setup([protein, new_mutation] * ":", "Hsa/jun2011")
    assert !mutation.truncated



  end


end


