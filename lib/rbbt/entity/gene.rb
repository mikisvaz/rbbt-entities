require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/sources/entrez'
require 'rbbt/sources/matador'
require 'rbbt/sources/cancer'
require 'rbbt/entity/protein'
require 'rbbt/entity/pmid'
require 'rbbt/entity/transcript'
require 'rbbt/bow/bow'

Workflow.require_workflow "Translation"

module Gene
  extend Entity

  def self.ensg2enst(organism, gene)
    @@ensg2enst ||= {}
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

  def self.gene_list_bases(genes)
    genes = genes.ensembl
    chromosome_genes = {}
    Misc.process_to_hash(genes){|genes| genes.chromosome}.each{|gene, chr| chromosome_genes[chr] ||= []; chromosome_genes[chr] << gene}
    total = 0
    chromosome_genes.each do |chr,gs|
      next if chr.nil?
      total += Misc.total_length(genes.annotate(gs).chr_range.compact)
    end
    
    total
  end

  def self.gene_list_exon_bases(genes)
    genes = genes.ensembl
    chromosome_genes = {}
    Misc.process_to_hash(genes){|genes| genes.chromosome}.each{|gene, chr| chromosome_genes[chr] ||= []; chromosome_genes[chr] << gene}

    @@exon_range_tsv ||= {}
    organism = genes.organism
    @@exon_range_tsv[organism] ||= Organism.exons(organism).tsv :persist => true, :fields => ["Exon Chr Start", "Exon Chr End"], :type => :list, :cast => :to_i
    total = 0

    chromosome_genes.each do |chr,gs|
      next if chr.nil?
      exons = genes.annotate(gs).transcripts.compact.flatten.exons.compact.flatten.uniq

      exon_ranges = exons.collect{|exon|
        Log.low "Exon #{ exon } does not have range" unless @@exon_range_tsv[organism].include? exon
        next unless @@exon_range_tsv[organism].include? exon
        pos = @@exon_range_tsv[organism][exon]
        (pos.first..pos.last)
      }.compact
      total += Misc.total_length(exon_ranges)
    end
    
    total
  end



  self.annotation :format
  self.annotation :organism

  self.format = Organism.identifiers("Hsa").all_fields - ["Ensembl Protein ID", "Ensembl Transcript ID"]

  property :ortholog => :array2single do |other|
    return self if organism =~ /^#{ other }(?!\w)/
    new_organism = organism.split(":")
    new_organism[0] = other
    new_organism = new_organism * "/"
    Gene.setup(Organism[organism]["ortholog_#{other}"].tsv(:persist => true).values_at(*self.ensembl).collect{|l| l.first}, "Ensembl Gene ID", new_organism)
  end
  persist :ortholog 

  property :to => :array2single do |new_format|
    return self if format == new_format
    genes = Translation.job(:tsv_translate, "", :organism => organism, :genes => self, :format => new_format).exec.values_at(*self)
    Gene.setup(genes, new_format, organism)
    genes
  end

  property :strand => :array2single do 
    @@strand_tsv ||= {}
    @@strand_tsv[organism] ||= Organism.gene_positions(organism).tsv(:fields => ["Strand"], :type => :single, :persist => true)
    to("Ensembl Gene ID").collect do |gene|
      @@strand_tsv[organism][gene]
    end
  end
  persist :_ary_strand

  property :ensembl => :array2single do
    to "Ensembl Gene ID"
  end

  property :biotype => :array2single do
    Organism.gene_biotype(organism).tsv(:persist => true, :type => :single).values_at *self.ensembl
  end
  persist :biotype

  property :entrez => :array2single do
    to "Entrez Gene ID"
  end

  property :uniprot => :array2single do
    to "UniProt/SwissProt Accession"
  end

  property :name => :array2single do
    return self if self.format == "Associated Gene Name"
    to "Associated Gene Name"
  end

  property :chr_start => :array2single do
    Organism.gene_positions(organism).tsv(:persist => true, :type => :single, :cast => :to_i, :fields => ["Gene Start"]).values_at *self
  end
  persist :chr_start

  property :go_bp_terms => :array2single do
    Organism.gene_go_bp(organism).tsv(:persist => true, :key_field => "Ensembl Gene ID", :fields => ["GO ID"], :type => :flat).values_at *self.ensembl
  end
  persist :go_bp_terms

  property :long_name => :array2single do
    entre = self.entrez
    gene = Entrez.get_gene(entrez).values_at(*entrez).collect{|gene| gene.nil? ? nil : gene.description.flatten.first}
  end
  persist :long_name

  property :description => :single2array do
    gene = Entrez.get_gene(to("Entrez Gene ID"))
    gene.nil? ? nil : gene.summary.flatten.first
  end
  persist :description

  property :transcripts => :array2single do
    res = Gene.ensg2enst(organism, self.ensembl)
    Transcript.setup(res, "Ensembl Transcript ID", organism)
    res
  end
  persist :transcripts

  property :proteins  => :array2single do
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
  persist :proteins

  property :max_transcript_length => :array2single do
    transcripts.collect{|list| list.sequence_length.compact.max}
  end
  persist :max_transcript_length

  property :max_protein_length => :array2single do
    proteins = self.proteins
    all_proteins = Protein.setup(proteins.flatten, "Ensembl Protein ID", organism)
    lengths = Misc.process_to_hash(all_proteins){|list| list.sequence_length}
    proteins.collect{|list| lengths.values_at(*list).compact.max}
  end
  persist :max_protein_length

  property :chromosome => :array2single do
    @@chromosome_tsv ||= {}
    @@chromosome_tsv[organism] ||= Organism.gene_positions(organism).tsv :fields => ["Chromosome Name"], :type => :single, :persist => true, :unnamed => true
    if Array === self
      to("Ensembl Gene ID").collect do |gene|
        @@chromosome_tsv[organism][gene]
      end
    else
      @@chromosome_tsv[organism][to("Ensembl Gene ID")]
    end
  end
  persist :chromosome

  property :chr_range => :array2single do
    chr_range_index ||= Organism.gene_positions(organism).tsv :fields => ["Gene Start", "Gene End"], :type => :list, :persist => true, :cast => :to_i
    to("Ensembl Gene ID").collect do |gene|
      next if not chr_range_index.include? gene
      Range.new *chr_range_index[gene]
    end
  end
  persist :chr_range

  property :articles => :array2single do
    PMID.setup(Organism.gene_pmids(organism).tsv(:persist => true, :fields => ["PMID"], :type => :flat, :unnamed => true).values_at *self.entrez)
  end
  persist :articles

  property :sequence => :array2single do
    @@sequence_tsv ||= {}
    @@sequence_tsv[organism] ||= Organism.gene_sequence(organism).tsv :persist => true, :unnamed => true
    @@sequence_tsv[organism].values_at *self.ensembl
  end
  persist :sequence

  property :matador_drugs => :array2single do
    @@matador ||= Matador.protein_drug.tsv(:persist => false).tap{|o| o.unnamed = true}

    ensg = self.to("Ensembl Gene ID")

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
  persist :matador_drugs

  property :drugs => :array2single do
    @matador_drugs = matador_drugs
  end
  persist :drugs

  property :kegg_pathway_drugs => :array2single do
    self.collect{|gene|
      pathway_genes = gene.kegg_pathways
      next if pathway_genes.nil?
      pathway_genes = pathway_genes.compact.flatten.genes.flatten
      Gene.setup(pathway_genes, "KEGG Gene ID", organism)

      pathway_genes.compact.drugs.compact.flatten.uniq
    }
  end
  persist :kegg_pathway_drugs

  property :pathway_drugs => :array2single do
    kegg_pathway_drugs
  end
  persist :pathway_drugs

  property :related_cancers => :array2single do
    Cancer["cancer_genes.tsv"].tsv(:persist => true, :type => :list).values_at(*self.name).collect{|v| v.nil? ? nil : (v["Tumour Types (Somatic Mutations)"].split(", ") + v["Tumour Types (Germline Mutations)"].split(", ")).uniq}
  end
  persist :related_cancers

  property :somatic_snvs => :array2single do
    names = self.name
    raise "No organism defined" if self.organism.nil?
    clean_organism = self.organism.sub(/\/.*/,'') + '/jun2011'
    names.organism = clean_organism
    ranges = names.chromosome.zip(name.chr_range).collect do |chromosome, range|
      next if range.nil?
      [chromosome, range.begin, range.end] * ":"
    end
    Sequence.job(:somatic_snvs_at_genomic_ranges, File.join("Gene", (names.compact.sort * ", ")[0..80]), :organism => clean_organism, :ranges  => ranges).fork.join.load.values_at *ranges
  end
  persist :somatic_snvs


  property :literature_score do |terms|
    terms = terms.collect{|t| t.stem}
    articles = self.articles
    if articles.nil? or articles.empty?
      0
    else
      articles.inject(0){|acc,article| acc += article.text.words.select{|word| terms.include? word}.length }.to_f / articles.length
    end
  end
  persist :literature_score


  property :ihop_interactions => :single do
    uniprot = self.uniprot
    if uniprot.nil?
      nil
    else
      sentences = []

      begin
        url = "http://ws.bioinfo.cnio.es/iHOP/cgi-bin/getSymbolInteractions?ncbiTaxId=9606&reference=#{uniprot}&namespace=UNIPROT__AC" 
        doc = Nokogiri::XML(Open.read(url))
        sentences = doc.css("iHOPsentence")
      rescue
      end

      sentences
    end
  end

  property :tagged_ihop_interactions => :single do
    interactors = []
    ihop_interactions = self.ihop_interactions
    if ihop_interactions.nil?
      nil
    else
      ihop_interactions.each do |sentence|
        sentence.css('iHOPatom').collect{|atom|
          atom.css('evidence');
        }.compact.flatten.each do |evidence|
          symbol =  evidence.attr('symbol')
          taxid  =  evidence.attr('ncbiTaxId')

          if Organism.entrez_taxids(self.organism).list.include? taxid
            interactors << symbol
          end
        end
      end

      Gene.setup(interactors, "Associated Gene Name", self.organism).organism
      
      interactors_ensembl = interactors.ensembl

      interactors2ensembl = {}
      interactors.collect{|i| i}.zip(interactors_ensembl.collect{|i| i}).each do |o,e|
        interactors2ensembl[o] = e
      end

      ihop_interactions.collect do |sentence|
        sentence.css('iHOPatom').each{|atom|
          literal = atom.content()
          evidences = atom.css('evidence')
          symbol = evidences.collect do |evidence|
            symbol =  evidence.attr('symbol')
            taxid  =  evidence.attr('ncbiTaxId')

            if Organism.entrez_taxids(self.organism).list.include? taxid
              symbol
            else
              nil
            end
          end.compact.first

          evidences.remove

          if interactors2ensembl.include? symbol and not interactors2ensembl[symbol].nil?
            atom.children.remove
            interactor = interactors2ensembl[symbol]
            atom.replace interactor.respond_to?(:link)? interactor.link(nil, nil, :html_link_extra_attrs => "title='#{literal}'") : interactor.name
          end
        }
        sentence.to_s
      end
    end
  end
end
