from doepy import build
from PyQt5 import QtWidgets

class DoEbuilder(QtWidgets.QWidget):
    def __init__(self, parent, main):
        super(DoEbuilder, self).__init__(parent)

        self.main = main

    def buildCC(self, bounds, variables):

        DoE_dict = {}
        for variable, bound in zip(variables, bounds):
            # new_dict_item = {variable : bound}
            DoE_dict[variable] = bound
            # DoE_dict = {**DoE_dict, **new_dict_item}

        design = build.central_composite(
            DoE_dict,
            center=(1, 1),
            alpha='orthogonal',
            face='cci'
        )

        return design