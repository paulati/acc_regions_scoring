from Bio import AlignIO
from os import listdir
import numpy as np
from maf_consensus import MafConsensus
from config import MafConsensusConfig, MafEncodeScoreConfig
from os.path import isfile, join
from scipy.spatial import distance


class MafEncodeScore:

    def __init__(self, ingroup_maf_base_path, chromosome, shift_specie):
        only_files = [f for f in listdir(ingroup_maf_base_path) if isfile(join(ingroup_maf_base_path, f))]
        print("cantidad de archivos maf a procesar:")
        print(str(len(only_files)))
        self.maf_files = {}
        for maf_file_name in only_files:
            self.maf_files[maf_file_name] = join(ingroup_maf_base_path, maf_file_name)
        self.__load_consensus(chromosome)
        self.shift_specie = shift_specie

    def __load_consensus(self, chromosome):
        maf_consensus = MafConsensus(chromosome)

        maf_consensus_config = MafConsensusConfig()

        # load_consensus_dict(self, consensus_base_path, consensus_file_name_head, consensus_file_name_tail)
        outgroup_consensus_base_path = maf_consensus_config.outgroup_consensus_base_path()
        outgroup_consensus_file_name_head = maf_consensus_config.outgroup_consensus_file_name_head()
        outgroup_consensus_file_name_tail = maf_consensus_config.outgroup_consensus_file_name_tail()
        self.outgroup_consensus_dict = maf_consensus.load_consensus_dict(outgroup_consensus_base_path,
                                                                         outgroup_consensus_file_name_head,
                                                                         outgroup_consensus_file_name_tail)

        ingroup_consensus_base_path = maf_consensus_config.ingroup_consensus_base_path()
        ingroup_consensus_file_name_head = maf_consensus_config.ingroup_consensus_file_name_head()
        ingroup_consensus_file_name_tail = maf_consensus_config.ingroup_consensus_file_name_tail()
        self.ingroup_consensus_dict = maf_consensus.load_consensus_dict(ingroup_consensus_base_path,
                                                                        ingroup_consensus_file_name_head,
                                                                        ingroup_consensus_file_name_tail)

    def process(self, out_file_path):

        result = {}

        file = open(out_file_path, "w")

        print("cantidad de archivos a procesar:")
        print(len(self.maf_files.keys()))

        for maf_file_name in self.maf_files.keys():

            maf_file_path = self.maf_files[maf_file_name]

            # TODO: cambiar este return a struct o clase
            (cols_sc, block_sc, shifts_c, hamming_in, hamming_out, hamming_ratio) = self.__process_maf_file(
                maf_file_name)

            result[maf_file_path] = (cols_sc, block_sc, shifts_c, hamming_in, hamming_out, hamming_ratio)

        for k in result.keys():
            line = k + "\t" + str(result[k][0]) + "\t" + str(result[k][1]) + "\t" + str(result[k][2]) + \
                   "\t" + str(result[k][3]) + "\t" + str(result[k][4]) + "\t" + str(result[k][5]) + "\n"
            #line = k + "\t" + str(result[k][2]) + \
            #      "\t" + str(result[k][3]) + "\t" + str(result[k][4]) + "\t" + str(result[k][5]) + "\n"
            # print(line)
            file.write(line)

        file.close()

    # def __count_shifts(self, alignment, maf_file_name):
    @staticmethod
    def __count_shifts(outgroup_consensus_seq, ingroup_consensus_seq, shift_seq):

        shifts = MafEncodeScore.__get_shifts(shift_seq, outgroup_consensus_seq)

        ingroup_shifts = MafEncodeScore.__get_shifts(shift_seq, ingroup_consensus_seq)

        result_shifts = np.zeros(len(shifts))

        for i in range(len(shifts)):
            elem_shift = shifts[i]
            elem_consensus_ingroup_shift = ingroup_shifts[i]
            # print(elem_shift)
            # print(elem_consensus_ingroup_shift)

            result_shifts[i] = elem_shift == 1 and elem_consensus_ingroup_shift == 0
            #if result_shifts[i] == 1:
            #    print("uno que sirve")

        shifts_count = np.sum(result_shifts)

        # print(shifts)
        # print(shifts_count)

        return shifts_count

    @staticmethod
    def __get_shifts(shift_sequence, consensus_sequence):
        result = None
        if len(shift_sequence) == len(consensus_sequence):
            result = np.zeros(len(shift_sequence))
            # print(shift_sequence)
            # print(consensus_sequence)
            for i in range(len(shift_sequence)):
                base_shift_sequence = shift_sequence[i]
                base_consensus_sequence = consensus_sequence[i]
                # TODO: sacar el hardcode de la X y pasarla a configuracion
                if base_consensus_sequence == "X":
                    result[i] = 0
                else:
                    if base_consensus_sequence == base_shift_sequence:
                        result[i] = 0
                    else:
                        result[i] = 1
        else:
            print("error, las secuencias deben ser de la misma longitud")
        return result

    @staticmethod
    def __get_specie_sequence(specie_name, alignment):

        result = ""

        for sequence in alignment:

            sequence_specie = sequence.id.split(".")[0]
            # print(sequence_specie)

            if sequence_specie == specie_name:
                result = sequence.seq

        return result

    @staticmethod
    def __score_columns(align_bit_mat):

        column_score = np.zeros(align_bit_mat.shape[1])

        for column_index in range(align_bit_mat.shape[1]):

            # print(align_bit_mat[:, column_index])

            # para cada columna, veo cuales son las filas que tienen un 0
            column = align_bit_mat[:, column_index]
            zero_indexes = np.where(column == 0)[0]
            # print(zero_indexes)

            # score column:

            col_score = 0
            from_index = 0
            for zero_index in zero_indexes:
                to_index = zero_index + 1
                # print("from: " + str(from_index))
                # print("to: " + str(to_index))
                block = column[from_index:to_index]  # incluye al limite inferior, NO incluye al limite superior
                block_sum = np.sum(block)
                block_len = len(block)

                # print(block)
                # print("len: " + str(block_len))
                # print("sum: " + str(block_sum))
                from_index = to_index
                col_score = col_score + float(block_sum) / float(block_len)
                # print("col_score: " + str(col_score))

            column_score[column_index] = col_score

        # print(column_score)

        return column_score

    @staticmethod
    def __score_block(columns_score):
        # score_block = np.sum(columns_score) / align_bit_mat.shape[1]
        score_block = np.sum(columns_score) / len(columns_score)
        # print(np.sum(columns_score))
        # print(len(columns_score))
        # print(score_block)
        return score_block

    @staticmethod
    def __alignment_bit_matrix(alignment_matrix):
        result = np.zeros(alignment_matrix.shape, dtype=int, order='C')
        row_count = alignment_matrix.shape[0]
        col_count = alignment_matrix.shape[1]
        for col_index in range(col_count):
            # la ultima fila no la miro y reviso hasta la de indice 0
            for row_index in range(row_count-2, -1, -1):
                # esto NO debe ser case sensitive:
                if(str.lower(str(alignment_matrix[row_index][col_index])) ==
                        str.lower(str(alignment_matrix[row_index+1][col_index]))):
                    result[row_index][col_index] = 1
                else:
                    result[row_index][col_index] = 0
        # print(result)
        return result

    @staticmethod
    def __alignment_matrix(alignment):

        rows_count = len(alignment)
        if rows_count > 0:
            cols_count = len(alignment[0].seq)
        else:
            cols_count = 0

        result = np.chararray((rows_count, cols_count))
        result[:] = '.'

        row_index = 0
        for seq_rec in alignment:
            result[row_index] = seq_rec.seq
            row_index = row_index + 1

        return result

    def __process_maf_file(self, maf_file_name):

        # print(maf_file_path)
        maf_file_path = self.maf_files[maf_file_name]

        alignment = AlignIO.read(maf_file_path, "maf")  # my input alignment
        align_matrix = MafEncodeScore.__alignment_matrix(alignment)
        align_bit_mat = MafEncodeScore.__alignment_bit_matrix(align_matrix)

        columns_score = MafEncodeScore.__score_columns(align_bit_mat)
        # print("columns_score: ")
        # print(columns_score)

        block_score = MafEncodeScore.__score_block(columns_score)
        # print("block_score: " + str(block_score))

        if maf_file_name in self.outgroup_consensus_dict.keys():
            outgroup_consensus_seq = self.outgroup_consensus_dict[maf_file_name]
        else:
            outgroup_consensus_seq = ""

        # print(outgroup_consensus_seq)

        if maf_file_name in self.ingroup_consensus_dict.keys():
            ingroup_consensus_seq = self.ingroup_consensus_dict[maf_file_name]
        else:
            ingroup_consensus_seq = ""

        # print(ingroup_consensus_seq)

        shift_seq = MafEncodeScore.__get_specie_sequence(self.shift_specie, alignment)
        # print(shift_seq)

        if len(outgroup_consensus_seq) > 0 and len(ingroup_consensus_seq) > 0 and len(shift_seq) > 0:

            # shift_specie = "ornAna1"
            shift_count = self.__count_shifts(outgroup_consensus_seq, ingroup_consensus_seq, shift_seq)
            # print("shift_count: " + str(shift_count))

            # print(ingroup_consensus_seq)
            # print(shift_seq)
            # print(outgroup_consensus_seq)

            ingroup_consensus_vec = list(ingroup_consensus_seq)
            outgroup_consensus_vec = list(outgroup_consensus_seq)

            hamming_ingroup = distance.hamming(ingroup_consensus_vec, shift_seq)

            hamming_outgroup = distance.hamming(outgroup_consensus_vec, shift_seq)

        else:
            shift_count = 0.0
            hamming_ingroup = 0.0
            hamming_outgroup = 0.0


        if hamming_ingroup == 0:
            if hamming_outgroup == 0:
                hamming_ratio = 0
            else:
                hamming_ratio = "MAX"
        else:
            hamming_ratio = hamming_outgroup / hamming_ingroup

        # jaccard_ingroup = distance.jaccard(ingroup_consensus_vec, shift_seq)
        # jaccard_outgroup = distance.jaccard(outgroup_consensus_vec, shift_seq)

        # TODO: cambiar este return a struct o clase
        return columns_score, block_score, shift_count, hamming_ingroup, hamming_outgroup, hamming_ratio


def main():

    chromosome = 22

    maf_encode_score_config = MafEncodeScoreConfig()

    ingroup_maf_base_path = maf_encode_score_config.ingroup_maf_base_path() + str(chromosome)

    shift_specie = maf_encode_score_config.shift_specie()

    maf_encode_score = MafEncodeScore(ingroup_maf_base_path, chromosome, shift_specie)

    out_file_path = maf_encode_score_config.out_file_path_head() + str(chromosome) + maf_encode_score_config.\
        out_file_path_tail()

    maf_encode_score.process(out_file_path)


# main()
