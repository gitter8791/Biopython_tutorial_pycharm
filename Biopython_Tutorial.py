'''
This is Biopython tutorial. PDF in project folder. 
'''

from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import GC

def f_location_test():
	my_snp = 4350
	record = SeqIO.read("NC_005816.gb", "genbank")
	for feature in record.features:
		if my_snp in feature:
			print("%s %s" % (feature.type, feature.qualifiers.get("db_xref")))


def f_seq_position():
	start_pos = SeqFeature.AfterPosition(5)
	end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
	my_location = SeqFeature.FeatureLocation(start_pos, end_pos)
	print(my_location, property(my_location), type(my_location))

	exact_location = SeqFeature.FeatureLocation(5, 9)
	exact_location.nofuzzy_start #??
	print(exact_location)


def f_record_objects_fasta():
	# record = SeqIO.read("NC_005816.fna", "fasta")
	record = SeqIO.read("NC_005816.gb", "genbank")

	print(record.id, '\n', record.name, '\n', record.description)
	print(repr(record.seq))
	print(len(record.annotations))
	print(record.annotations["source"])
	print(record.dbxrefs)
	print(len(record.features))


def f_create_sequence_record():
	simple_seq = Seq("GATC", )
	simple_seq_r = SeqRecord(simple_seq, id="AC12345")
	simple_seq_r.id = "AC12345"  # reduntant
	simple_seq_r.description = "Test seq for creating record"
	simple_seq_r.annotations["a_custom_field"] = "This is is annotation, string: a_custom_field"
	simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]  # not supported in BioSQL

	print(simple_seq_r.letter_annotations)
	# print(simple_seq_r.annotations)
	# print(simple_seq_r.annotations["a_custom_field"])
	print('simple_seq: ', simple_seq, '\n', simple_seq_r)


def f_seq_list_join():
	list_of_seqs = [Seq("ACGT", generic_dna), Seq("AACC", generic_dna), Seq("GGTT", generic_dna)]
	concatenated = Seq("", generic_dna)
	for s in list_of_seqs:
		concatenated += s

	print(concatenated)


# modules for testing Sequence handling functions

def f_basic_alphbet():
	my_seq = Seq("AGTACACTGGT")

	print('sequence:', my_seq)
	# AGTACACTGGT
	print('compementary: ', my_seq.complement())
	# TCATGTGACCA
	print('complementary reversed :', my_seq.reverse_complement())
	# ACCAGTGTACT
	print('type general: ', type(my_seq))
	# <class 'Bio.Seq.Seq'>
	print('general alphbet: ', my_seq.alphabet, '\n')
	# Alphabet()


# sequence handled as protein aminoacids
def f_seq_as_prot():
	my_prot = Seq("AGTACACTGGT", IUPAC.protein)
	print('protein sequence: ', my_prot)
	print('type protein seq: ', type(my_prot))
	print('protein alphabet: ', my_prot.alphabet)


# sequence handled as dna
def f_seq_hanadled_as_DNA():
	my_seq_IUPAC_DNA = Seq("AGTACACTGGTACGT", IUPAC.unambiguous_dna)
	print('dna sequence: ', my_seq_IUPAC_DNA)

	# sequences can be handled as strings
	print('dna sequence type: ', type(my_seq_IUPAC_DNA))
	print('dna alphbet: ', my_seq_IUPAC_DNA.alphabet)
	print('sequence length: ', len(my_seq_IUPAC_DNA))
	print('GC ratio in sequence: ',
	      100 * float(my_seq_IUPAC_DNA.count("G") + my_seq_IUPAC_DNA.count("C")) / len(my_seq_IUPAC_DNA))

	# built in function to determine GC ratio
	# from Bio.SeqUtils import GC
	print('GC ratio with b.in func: ', GC(my_seq_IUPAC_DNA))

	# print out seq items
	print('sequence items from 5 to 12: ', my_seq_IUPAC_DNA[4:12])
	print('first letter in every codon: ', my_seq_IUPAC_DNA[0::3])
	print('     sec letter: ', my_seq_IUPAC_DNA[1::3])
	print('     third letter: ', my_seq_IUPAC_DNA[2::3])

	print('back to string: ', type(str(my_seq_IUPAC_DNA)))  # str() redundant??

	print('\n')


	# for index, letter in enumerate(my_seq_IUPAC_DNA):
	#     print("index: %i; letter: %s" % (index, letter))
	# print('\n')
	#
	# print('first letter: ', my_seq_IUPAC_DNA[0])  # first letter
	# print('third letter: ', my_seq_IUPAC_DNA[2])  # third letter
	# print('last letter: ', my_seq_IUPAC_DNA[-1])  # last letter
	# print('\n')
	#
	# #non overlapping count
	# print('non-overlapping count (AA in AAA): ', Seq("AAAA").count("AA"))
	# #overlapping count??
	#


# takes fasta file and prints sequence data and id and length
def f_FASTA_tut():
	for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
		print(seq_record.id)
		print(repr(seq_record.seq))
		print(len(seq_record))


# takes gbk file and prints id, seq, and length
def f_GBank_file():
	for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
		print(seq_record.id)
		print(repr(seq_record.seq))  # repr() cuts sequence print out
		print(len(seq_record))


# main function. uncomment function calls for testing
def main():
	pass


	# f_record_objects_fasta()
	# f_create_sequence_record
	# f_basic_alphbet
	# f_GBank_file()
	# f_FASTA_tut()
	# f_seq_as_prot()
	# f_seq_list_join()
	# f_seq_hanadled_as_DNA()



if __name__ == '__main__':
	main()
