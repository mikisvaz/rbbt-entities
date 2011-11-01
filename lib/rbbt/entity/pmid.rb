require 'rbbt/entity'
require 'rbbt/sources/pubmed'

module PMID
  extend Entity

  self.format = "PMID"

  property :title => :array2single do
    @title ||= begin
                 PubMed.get_article(self).values_at(*self).collect{|article| article.nil? ? nil : article.title}
               end
  end

  property :pubmed_url => :single2array do
    "<a class='pmid' href='http://www.ncbi.nlm.nih.gov/pubmed/#{self}'>#{ self }</a>"
  end
end

