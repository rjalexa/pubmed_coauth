"""Write a GEXF graph file from coauthrships of PubMed."""
# Build a multigraph of coauthorships
# Launch a pubmed search on a given subject or author (list)
# Each node is the name of an author
# Each edge is a coauthorship between each couple of authors
# Each edge will have the following attributes
# from the BioPython medline interface
# * DP: date of publication YYYY (MM) (DD)
# * JT: Journal -> name of journal
# * PT: Paper type -> a list in which you can find 'journal article'
# * PL: Publication location -> country of location
# * PMID: Pubmed publication unique ID
# * AD: affiliation --- not used yet

from Bio import Entrez  # to access the pubmed database
from Bio import Medline
import itertools  # to generate all coauthors combinations
import networkx as nx  # to create, analyze and write the graph to a gexf format
import time  # to generate a time stamped gexf filename
from affiliation_parser import parse_affil  # depends on nltk and sklearn

# use a file to save results if file does not exist;
# if the file exists use it instead of fetching the data from pubmed;
# datafilename = "./author_correlation_file.txt"

# if os.path.isfile(datafilename):
#    records = []
#    with open(datafilename, 'r') as f:
#        for x in f:
#            records.extend(f.read(str(x)))
# else:   # there is no data saved in a file hence fetch it from PubMed
# setup your email for PubMed searches
Entrez.email = "rja@it.ibm.com"
# setup a couple of other pubmed search parameters
maxrecords = "100000"
corpus = "pubmed"

# the following are sample searchstring
# could also limit by year of publication (since 2014 authors have affiliations)
# searchstring='("lung diseases, interstitial"[MeSH Terms] NOT "pneumoconiosis"[MeSH Terms]) AND "artificial intelligence"[MeSH Terms] AND "humans"[MeSH Terms]'
# a toy search
# searchstring='alexander rj[Author]'
# a medium sized one
# searchstring='("Diabetes Mellitus"[Mesh]) AND "Artificial Intelligence"[Mesh]'
# another medium one
# searchstring='((("Pancreatitis"[Mesh] NOT "Pancreatitis, Alcoholic"[Mesh]) NOT "Pancreatitis, Chronic"[Mesh]) NOT "Pancreatitis, Graft"[Mesh]) AND Sphincterotomy, Endoscopic[Mesh] '
# searchstring = '(Santoro, Armando[Full Author Name] OR Santoro, Armando[Full Investigator Name]) AND (Scorsetti, Marta[Full Author Name] OR Scorsetti, Marta[Full Investigator Name])'
# searchstring = '(Scorsetti, Marta[Full Author Name] OR Scorsetti, Marta[Full Investigator Name])'
# searchstring = 'mastroianni cm[Author] OR Lichtner, Miriam[Full Author Name]'
searchstring = '"Artificial Intelligence"[Mesh]'


handle = Entrez.esearch(db=corpus, term=searchstring,
                        retmax=maxrecords, sort="relevance", retmode="xml")
records = Entrez.read(handle)

ids = records['IdList']
print("Found {} articles using the following query string: {}".format(
       records["Count"], searchstring))

h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')

# records is a generator and therefore can only be iterated once
records = Medline.parse(h)

#    # write all retrieved records to a file
#    with open(datafilename, 'w') as f: #using open within a with block automatically closes the file at exit
#        for x in records:
#            f.write(str(x))


# now that you have the records either from file or efetched go ahead
# create the new multigraph. A multigraph can hold N edges between a couple of nodes. It limits some algos
# but in this case I will have one edge for each coauthorship
# we will use authors_dict for the nodes, coa_df for the edges, pub_df for the edges attributes
G = nx.MultiGraph(name=searchstring)

# now records hold a dictionary of matching publications
# warning: if already used to write the file above it will be exhausted already
for record in records:
    au = record.get('AU', '?')  # list of authors
    dp = record.get('DP', '?')  # date of publication
    # keep only the first word in the string, the year
    yp = dp.split(' ', 1)[0]
    pi = record.get('PMID', '?')  # publication ID
    af = record.get('AD', '?')  # affiliation string
    affiliation = parse_affil(af)
    country = affiliation['country']  # country of affiliation of first author; can be blank;
    # loop through each paper's author list and add edges to the graph with auth1-auth2 and PMID as attribute
    coauth_pairs_ls = list(itertools.combinations(au, 2))
    for pair in coauth_pairs_ls:
        auth1, auth2 = pair
        G.add_edge(auth1, auth2, pmid=pi, year=yp, auth1_country=country)


# print some Graph statistics
print("Some statistics of this graph:")
print(nx.info(G))
print("The density of the graph is {}.".format(nx.density(G)))
print("The degree of the 'Alexander RJ' node is {}.".format(
    nx.degree(G, "Alexander RJ")))
# print(nx.degree_centrality(G))
# print(nx.betweenness_centrality(G))
# print(nx.eigenvector_centrality(G))
gclique = list(nx.find_cliques(G))
print("There are {} cliques in the graph.".format(len(gclique)))
# Find connected components
gcomps = nx.connected_components(G)
# gcomps is a generator so I am counting its output in a loop
numwords = sum(1 for w in gcomps)
print("There are {} connected component in the graph.".format(numwords))
# now save the graph to a GEXF file to use with Gephi
timestr = time.strftime("%Y%m%d-%H%M%S")
filename = "coauthors" + timestr + ".gexf"
print("Saving the graph in GEXF format to {} file.".format(filename))
nx.write_gexf(G, filename)
