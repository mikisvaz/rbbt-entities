require 'rbbt/entity'
require 'rbbt/sources/pubmed'

module PMID
  extend Entity

  self.format = "PMID"

  property :article => :array2single do
    @article ||= begin
                   PubMed.get_article(self).values_at(*self)
                 end
  end

  property :title => :array2single do
    @title ||= begin
                 article.collect{|a| a.nil? ? nil : a.title}
               end
  end

  property :text => :array2single do
    @text ||= begin
                 article.collect{|a| a.nil? ? nil : a.text}
               end
  end


  property :pubmed_url => :single2array do
    "<a class='pmid' href='http://www.ncbi.nlm.nih.gov/pubmed/#{self}'>#{ self }</a>"
  end
end

