import csv
import gzip
import sqlite3

from contextlib import closing


def create_string_db(db_filename, protein_links_filename, protein_aliases_filename):
    with closing(sqlite3.connect(db_filename)) as db:
        db.text_factory = str
    
        insert_protein_ids(db, protein_links_filename)
        insert_protein_links(db, protein_links_filename)
        insert_protein_names(db, protein_aliases_filename)


###
#
# Internal protein IDs
#

def create_proteins_table(db):
    db.execute(
        """
        CREATE TABLE proteins(
            _id INTEGER PRIMARY KEY,
            protein_id VARCHAR)
        """)


def insert_protein_ids(db, filename):
    create_proteins_table(db)
    
    with gzip.open(filename) as stream:
        reader = csv.reader(stream, csv.excel_tab, delimiter=" ")

        header = reader.next()
        check_protein_links_header(header)

        protein_ids = {}
        next_protein_id = 0
        
        for record in reader:
            if record[PPI_PROTEIN1_COLUMN].startswith("9606.") and int(record[PPI_SCORE_COLUMN]) > 800:
                species1, protein1 = record[PPI_PROTEIN1_COLUMN].split(".", 1)
                species2, protein2 = record[PPI_PROTEIN2_COLUMN].split(".", 1)
                assert species1 == species2

                for protein in [protein1, protein2]:
                    if protein not in protein_ids:
                        protein_ids[protein] = next_protein_id
                        next_protein_id += 1

        with db:
            for protein, protein_id in sorted(protein_ids.iteritems(), key=lambda x: x[1]):
                db.execute("INSERT INTO proteins VALUES (?, ?)", (protein_id, protein))


def read_protein_ids(db):
    return dict(db.execute("SELECT protein_id, _id FROM proteins"))


####
#
# Protein links
#

PPI_PROTEIN1_COLUMN = 0
PPI_PROTEIN2_COLUMN = 1
PPI_SCORE_COLUMN = 9


def check_protein_links_header(header):
    assert header[PPI_PROTEIN1_COLUMN] == "protein1"
    assert header[PPI_PROTEIN2_COLUMN] == "protein2"
    assert header[PPI_SCORE_COLUMN] == "combined_score"


def create_protein_links_table(db):
    db.execute(
        """
        CREATE TABLE protein_links(
            _id INTEGER PRIMARY KEY,
            protein_id_a INTEGER REFERENCES proteins,
            protein_id_b INTEGER REFERENCES proteins,
            combined_score INTEGER)
        """)

    db.execute("CREATE INDEX protein_link_index on protein_links(protein_id_a)")


def insert_protein_links(db, filename):
    create_protein_links_table(db)
    
    with gzip.open(filename) as stream:
        reader = csv.reader(stream, csv.excel_tab, delimiter=" ")

        header = reader.next()
        check_protein_links_header(header)

        protein_ids = read_protein_ids(db)
        record_id = 0
        
        with db:
            for record in reader:
                if record[PPI_PROTEIN1_COLUMN].startswith("9606.") and int(record[PPI_SCORE_COLUMN]) > 800:
                    species1, protein1 = record[PPI_PROTEIN1_COLUMN].split(".", 1)
                    species2, protein2 = record[PPI_PROTEIN2_COLUMN].split(".", 1)
                    assert species1 == species2
                    score = int(record[PPI_SCORE_COLUMN])

                    db.execute(
                        "INSERT INTO protein_links VALUES(?, ?, ?, ?)",
                        (record_id, protein_ids[protein1], protein_ids[protein2], score))

                    record_id += 1


####
#
# Protein names
#

NAMES_PROTEIN_COLUMN = 0
NAMES_ALIAS_COLUMN = 1
NAMES_SOURCE_COLUMN = 2


def check_names_header(header):
    assert header[0] == "## string_protein_id ## alias ## source ##"


def create_names_table(db):
    db.execute(
        """
        CREATE TABLE protein_names(
            protein_id INTEGER REFERENCES proteins,
            protein_name VARCHAR,
            source VARCHAR,
            PRIMARY KEY (protein_id, protein_name, source))
        """)

    db.execute("CREATE INDEX protein_name_index ON protein_names(protein_name)")


def insert_protein_names(db, filename):
    create_names_table(db)
    
    with gzip.open(filename) as stream:
        reader = csv.reader(stream, csv.excel_tab)

        header = reader.next()
        check_names_header(header)

        protein_ids = read_protein_ids(db)
        
        with db:
            for record in reader:
                if record[NAMES_PROTEIN_COLUMN].startswith("9606."):
                    species, protein = record[NAMES_PROTEIN_COLUMN].split(".", 1)
                    alias = record[NAMES_ALIAS_COLUMN]
                    sources = record[NAMES_SOURCE_COLUMN].split()

                    if protein in protein_ids:
                        for source in sources:
                            if source == "BioMart_HUGO":
                                db.execute(
                                    "INSERT INTO protein_names VALUES (?, ?, ?)",
                                    (protein_ids[protein], alias, source))
