require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/sources/entrez'
require 'rbbt/sources/matador'
require 'rbbt/entity/protein'
require 'rbbt/entity/pmid'

Workflow.require_workflow "Translation"

module Gene
  extend Entity

  def self.ensg2enst(organism, gene)
    @@ensg2enst ||= {}
    #@@ensg2enst[organism] ||= Organism.gene_transcripts(organism).tsv(:type => :single, :key_field => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true).tap{|o| o.unnamed = true}
    @@ensg2enst[organism] ||= Organism.gene_transcripts(organism).tsv(:type => :flat, :key_field => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true).tap{|o| o.unnamed = true}

    if Array === gene
      @@ensg2enst[organism].values_at *gene
    else
      @@ensg2enst[organism][gene]
    end
  end


  def self.filter(query, field = nil, options = nil, entity = nil)
    return true if query == entity
    
    return true if query == Gene.setup(entity.dup, options.merge(:format => field)).name

    false
  end

  self.annotation :format
  self.annotation :organism

  self.format = Organism::Hsa.identifiers.all_fields - ["Ensembl Protein ID", "Ensembl Transcript ID"]

  property :ortholog => :array2single do |other|
    return self if organism =~ /^#{ other }(?!\w)/
    new_organism = organism.split(":")
    new_organism[0] = other
    new_organism = new_organism * "/"
    Gene.setup(Organism[organism]["ortholog_#{other}"].tsv(:persist => true).values_at(*self.ensembl).collect{|l| l.first}, "Ensembl Gene ID", new_organism)
  end

  property :to! => :array2single do |new_format|
    return self if format == new_format
    Gene.setup(Translation.job(:tsv_translate, "", :organism => organism, :genes => self, :format => new_format).exec.values_at(*self), new_format, organism)
  end

  property :to => :array2single do |new_format|
    return self if format == new_format
    to!(new_format).collect!{|v| Array === v ? v.first : v}
  end

  property :strand => :array2single do 
    @strand ||= Organism.gene_positions(organism).tsv(:fields => ["Strand"], :type => :single, :persist => true).values_at *self
  end

  property :ensembl => :array2single do
    @ensembl ||= to "Ensembl Gene ID"
  end

  property :entrez => :array2single do
    @entrez ||= to "Entrez Gene ID"
  end

  property :uniprot => :array2single do
    @uniprot ||= to "UniProt/SwissProt Accession"
  end

  property :name => :array2single do
    @name ||= to "Associated Gene Name"
  end

  property :chr_start => :array2single do
    @chr_start = begin
                   Organism.gene_positions(organism).tsv(:persist => true, :type => :single, :cast => :to_i, :fields => ["Gene Start"]).values_at *self
                 end
  end

  property :go_bp_terms => :array2single do
    @go_bp_terms ||= Organism.gene_go_bp(organism).tsv(:persist => true, :key_field => "Ensembl Gene ID", :fields => ["GO ID"], :type => :flat).values_at *self.ensembl
  end

  property :long_name => :single2array do
    gene = Entrez.get_gene(to("Entrez Gene ID"))
    gene.nil? ? nil : gene.description.flatten.first
  end

  property :description => :single2array do
    gene = Entrez.get_gene(to("Entrez Gene ID"))
    gene.nil? ? nil : gene.summary.flatten.first
  end

  property :transcripts => :array2single do
    res = Gene.ensg2enst(organism, self.ensembl)
    Transcript.setup(res, "Ensembl Transcript ID", organism)
    res
  end

  property :proteins  => :array2single do
    @proteins ||= begin
                    transcripts = Gene.ensg2enst(organism, self.ensembl)

                    all_transcripts = Transcript.setup(transcripts.flatten.compact.uniq, "Ensembl Transcript ID", organism)

                    transcript2protein = Misc.process_to_hash(all_transcripts){|list|
                      list.protein
                    }

                    res = transcripts.collect{|list|
                      Protein.setup(transcript2protein.values_at(*list).compact.uniq, "Ensembl Protein ID", organism)
                    }

                    Protein.setup(res, "Ensembl Protein ID", organism)
                  end
  end

  property :max_transcript_length => :array2single do
    transcripts.collect{|list| list.sequence_length.compact.max}
  end

  property :max_protein_length => :array2single do
    @max_protein_length ||= begin
                              proteins = self.proteins
                              all_proteins = Protein.setup(proteins.flatten, "Ensembl Protein ID", organism)
                              lengths = Misc.process_to_hash(all_proteins){|list| list.sequence_length}
                              proteins.collect{|list| lengths.values_at(*list).compact.max}
                            end
  end

  property :chromosome => :array2single do
    chr = Organism.gene_positions(organism).tsv :fields => ["Chromosome Name"], :type => :single, :persist => true
    chr.unnamed = true
    if Array === self
      to("Ensembl Gene ID").collect do |gene|
        chr[gene]
      end
    else
      chr[to("Ensembl Gene ID")]
    end
  end

  property :range => :array2single do
    pos = Organism.gene_positions(organism).tsv :fields => ["Gene Start", "Gene End"], :type => :list, :persist => true, :cast => :to_i
    to("Ensembl Gene ID").collect do |gene|
      next if not pos.include? gene
      Range.new *pos[gene]
    end
  end

  property :articles => :array2single do
    @articles ||= begin
                    PMID.setup(Organism.gene_pmids(organism).tsv(:persist => true, :fields => ["PMID"], :type => :flat).values_at *self.entrez)
                  end 
  end

  property :sequence => :array2single do
    @gene_sequence ||= Organism.gene_sequence(organism).tsv :persist => true
    @gene_sequence.unnamed = true
    @gene_sequence.values_at *self.ensembl
  end

  property :matador_drugs => :array2single do
    @matador_drugs ||= begin
                         @@matador ||= Matador.protein_drug.tsv(:persist => false).tap{|o| o.unnamed = true}
                         
                           ensg = self._to("Ensembl Gene ID")

                           transcripts = Gene.ensg2enst(organism, ensg)

                           t2ps = Misc.process_to_hash(transcripts.compact.flatten.uniq){|l| Transcript.enst2ensp(organism, l).flatten.compact.uniq}

                           all_proteins = t2ps.values.flatten.compact

                           chemical_pos = @@matador.identify_field "Chemical"

                           p2ds = Misc.process_to_hash(all_proteins){|proteins| 
                             @@matador.values_at(*proteins).collect{|values| 
                               next if values.nil?
                               values[chemical_pos]
                             }
                           }

                           res = transcripts.collect do |ts|
                             ps = t2ps.values_at(*ts).compact.flatten
                             p2ds.values_at(*ps).flatten.compact.uniq
                           end

                           res
                         end
  end

  property :drugs => :array2single do
    @matador_drugs = matador_drugs
  end

  property :kegg_pathway_drugs => :array2single do
    @kegg_patyhway_drugs ||= begin
                               self.collect{|gene|
                                 pathway_genes = gene.kegg_pathways
                                 next if pathway_genes.nil?
                                 pathway_genes = pathway_genes.compact.flatten.genes.flatten
                                 Gene.setup(pathway_genes, "KEGG Gene ID", organism)

                                 pathway_genes.compact.drugs.compact.flatten.uniq
                               }
                             end
  end

  property :pathway_drugs => :array2single do
    @keep_pathway_drugs ||= kegg_pathway_drugs
  end


end

module Transcript
  extend Entity

  self.annotation :format
  self.annotation :organism

  self.format = "Ensembl Transcript ID"

  def self.enst2ensg(organism, transcript)
    @@enst2ensg ||= {}
    @@enst2ensg[organism] ||= Organism.gene_transcripts(organism).tsv(:type => :single, :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Gene ID"], :persist => true).tap{|o| o.unnamed = true}
    res = if Array === transcript
            @@enst2ensg[organism].values_at *transcript
          else
            @@enst2ensg[organism][transcript]
          end
    Gene.setup(res, "Ensembl Gene ID", organism)
  end

  def self.enst2ensp(organism, transcript)
    @@enst2ensp ||= {}
    @@enst2ensp[organism] ||= Organism.transcripts(organism).tsv(:type => :single, :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Protein ID"], :persist => true)
    res = if Array === transcript
            @@enst2ensp[organism].values_at *transcript
          else
            @@enst2ensp[organism][transcript]
          end
    Protein.setup(res, "Ensembl Protein ID", organism)
  end


  property :to! => :array2single do |new_format|
    return self if format == new_format
    Gene.setup(Translation.job(:tsv_translate_probe, "", :organism => organism, :probes => self, :format => new_format).exec.values_at(*self), new_format, organism)
  end

  property :to => :array2single do |new_format|
    return self if format == new_format
    to!(new_format).collect!{|v| v.nil? ? nil : v.first}
  end

  property :ensembl => :array2single do
    to "Ensembl Transcript ID"
  end

  property :sequence => :array2single do
    transcript_sequence = Organism.transcript_sequence(organism).tsv :persist => true
    transcript_sequence.unnamed = true
    transcript_sequence.values_at *self.ensembl
  end

  property :sequence_length => :array2single do
    sequence.collect{|s|
      s.nil? ? nil : s.length
    }
  end

  property :gene => :array2single do
    Transcript.enst2ensg(organism, self)
  end

  property :protein => :array2single do
    Transcript.enst2ensp(organism, self)
  end
end

