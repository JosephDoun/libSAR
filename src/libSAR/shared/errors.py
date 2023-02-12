
class InvalidInputError(IOError):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
        
        