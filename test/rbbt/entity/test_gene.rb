require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')

require 'test/unit'
require 'rbbt/util/tmpfile'
require 'test/unit'
require 'rbbt/entity/gene'

class TestGene < Test::Unit::TestCase
  CDK5 = Gene.setup("CDK5", "Associated Gene Name", "Hsa")
  TP53 = Gene.setup("TP53", "Associated Gene Name", "Hsa")
  TWO = Gene.setup(["CDK5", "TP53"], "Associated Gene Name", "Hsa")

  def test_to
    assert_equal "1020", Gene.setup("CDK5", "Associated Gene Name", "Hsa").to("Entrez Gene ID")
  end

  def test_long_name
    assert_equal "cyclin-dependent kinase 5", Gene.setup("CDK5", "Associated Gene Name", "Hsa").long_name
    assert_equal ["cyclin-dependent kinase 5"], Gene.setup(["CDK5"], "Associated Gene Name", "Hsa").long_name

    assert_match /tumor/, Gene.setup("TP53", "Associated Gene Name", "Hsa").description
  end

  def test_transcripts
    assert CDK5.transcripts.length > 1
    assert_equal "Hsa", CDK5.transcripts.organism
    assert_equal "Hsa", CDK5.make_list.transcripts[0].organism
    assert_equal "Hsa", CDK5.transcripts.make_list.organism
  end

  def test_proteins
    assert CDK5.proteins.length > 1
  end

  def test_max_protein_length
    assert CDK5.max_protein_length > 200
    assert Array === TWO.max_protein_length
    assert TWO.max_protein_length.first > 200
  end

  def test_max_transcript_length
    assert CDK5.max_transcript_length > 200
    assert Array === TWO.max_transcript_length
    assert TWO.max_transcript_length.first > 200
  end

  def test_range
    assert Range === CDK5.chr_range
    assert Range === CDK5.make_list.chr_range.first
  end
end


