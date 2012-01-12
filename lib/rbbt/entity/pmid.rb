require 'rbbt/entity'
require 'rbbt/sources/pubmed'

module PMID
  extend Entity

  self.format = "PMID"

  property :article => :array2single do
    PubMed.get_article(self).values_at(*self)
  end
  persist :article

  property :title => :array2single do
    article.collect{|a| a.nil? ? nil : a.title}
  end
  persist :title

  property :text => :array2single do
    article.collect{|a| a.nil? ? nil : a.text}
  end
  persist :text 

  property :pubmed_url => :single2array do
    "<a class='pmid' href='http://www.ncbi.nlm.nih.gov/pubmed/#{self}'>#{ self }</a>"
  end
  persist :pubmed_url 
end

