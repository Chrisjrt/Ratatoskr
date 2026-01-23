import Entrez

def entrez_search(db, term, email, retmax=10000):
    with Entrez.esearch(db=db, term=term, email=email, retmax=retmax) as h:
        return Entrez.read(h)
    