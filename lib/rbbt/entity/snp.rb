require 'rbbt/entity'

module SNP
  extend Entity

  self.format = ["SNP", "SNP ID", "RSID"]

  def self.dbSNP_info
    @@dbSNP_info ||= DbSNP.mutations.tsv :persist => true
  end

  property :dbSNP_info => :array2single do
    SNP.dbSNP_info.chunked_values_at self
  end
end
