from converter.datamodel.submittable import AccessionedSubmittable


class Protocol(AccessionedSubmittable):

    def __init__(self, **kwargs):
        AccessionedSubmittable.__init__(self, **kwargs)
        self.name = kwargs.get("name")
        self.protocol_type = kwargs.get("protocol_type")
        self.hardware = kwargs.get("hardware")
        self.software = kwargs.get("software")

