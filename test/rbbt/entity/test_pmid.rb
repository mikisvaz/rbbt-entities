require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')

require 'test/unit'
require 'rbbt/util/tmpfile'
require 'test/unit'
require 'rbbt/entity/pmid'

class TestProtein < Test::Unit::TestCase
  
  def test_pmid_id
    assert_match /^PMID/, PMID.setup("21904853").docid
  end

  def test_pmid_text
    assert_match /TET2/, PMID.setup("21904853").text
  end
end


