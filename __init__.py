import DMCModel
import HeterogeneousGoModel
import HomogeneousGoModel

def get_model(type):
    if type == "HomGo":
        model = HomogeneousGoModel.HomogeneousGoModel()
    elif type == "HetGo":
        model = HeterogeneousGoModel.HeterogeneousGoModel()
    elif type == "DMC":
        model = DMCModel.DMCModel()
    else:
        print "ERROR. No such model exists."
    return model
