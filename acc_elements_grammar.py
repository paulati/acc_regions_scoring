from pyparsing import *


class ElemScoreGrammar:

    def __init__(self):

        self.acc_element_file_path_head = Literal("/home/rstudio/disco_tmp/data/acc_maf_ingroup/")

        self.acc_element_file_path = Combine(self.acc_element_file_path_head + OneOrMore(Word(printables))
                                             + WordEnd(".maf"))

        self.acc_element_score_num = Combine(Word(nums) + ZeroOrMore("." + Word(nums))) ^ \
                                     Combine(Word(nums) + Literal("."))
        self.acc_element_score_max = Literal("MAX")
        self.acc_element_score = self.acc_element_score_num ^ self.acc_element_score_max

        # acc_element_score_bases_matrix_row = Combine(OneOrMore(acc_element_score + matrix_separator))
        # acc_element_score_bases_matrix = Combine("[" + OneOrMore(acc_element_score_bases_matrix_row) + "]")

        self.acc_element_score_average = self.acc_element_score

        self.acc_element_shift_count = Combine(Word(nums) + "." + Word(nums))

        self.acc_element_hamming_ingroup = self.acc_element_score

        self.acc_element_hamming_outgroup = self.acc_element_score

        self.acc_element_hamming_relative = self.acc_element_score

        self.elem_last_line_tail = OneOrMore(self.acc_element_score) + Literal("]") + \
            self.acc_element_score_average.setResultsName("element_score_average") + \
            self.acc_element_shift_count.setResultsName("shift_count") + \
            self.acc_element_hamming_ingroup.setResultsName("hamming_ingroup") + \
            self.acc_element_hamming_outgroup.setResultsName("hamming_outgroup") + \
            self.acc_element_hamming_relative.setResultsName("hamming_relative")


class ElemScore:

    # TODO: checkear los tipos de datos de los atributos

    def __init__(self, elem_id, elem_last_line):

        grammar = ElemScoreGrammar()
        data = grammar .elem_last_line_tail.parseString(elem_last_line)
        self.elem_score_average = data["element_score_average"]
        self.shift_count = data["shift_count"]
        self.hamming_ingroup = data["hamming_ingroup"]
        self.hamming_outgroup = data["hamming_outgroup"]
        self.hamming_relative = data["hamming_relative"]
        self.elem_id = elem_id


class ElemScoreParser:

    def __init__(self, file_path):

        f = open(file_path, 'r')

        lines = f.readlines()
        lines_count = len(lines)

        element_line_index = {}

        grammar = ElemScoreGrammar()

        for index in range(lines_count):
            line = lines[index]
            elem_id = grammar.acc_element_file_path.searchString(line)
            if len(elem_id) > 0:
                element_line_index[elem_id] = index

        self.score_elements = []

        elems_count = len(element_line_index)

        for index in range(0, elems_count):

            elem_id = list(element_line_index.keys())[index]

            if index != (elems_count - 1):
                next_elem_line_index = list(element_line_index.values())[index + 1]
            else:
                next_elem_line_index = lines_count

            elem_last_line = lines[next_elem_line_index - 1]
            elem_score = ElemScore(elem_id, elem_last_line)
            self.score_elements.append(elem_score)

        f.close()
#
#
# filename = "/paula/2018/acc_regions/scoring/data/results/chr13_results_25.txt"
# elem_score_parser = ElemScoreParser(filename)
# result = elem_score_parser.score_elements
#
# print(result[0].elem_id)
# print(result[0].hamming_outgroup)
# print(result[0].shift_count)
# print(result[0].elem_score_average)
# print(result[0].hamming_ingroup)
# print(result[0].hamming_outgroup)
# print(result[0].hamming_relative)
#

