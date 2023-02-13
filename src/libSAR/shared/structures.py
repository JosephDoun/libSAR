from   abc    import ABC, abstractmethod
from   lxml   import etree
from   typing import List

import os


class SARImage(ABC):
    "Base class of SAR images defining abstract methods."
    
    _POLARISATIONS = {
            'S': {'H': ['HH'],
                  'V': ['VV']},
            "D": {'H': ['HH',
                        'HV'],
                  'V': ['VV',
                        'VH']}
            }
    
    def __init__(self, PATH: str) -> None:
        # super().__init__()
        if not os.path.exists(PATH):
            raise FileNotFoundError("File does not exist.")
    
    @property
    @abstractmethod
    def bands(self):
        ...


class XMLMetadata(ABC):
    "Base class for XML metadata. Base XML parser."
    __tag = lambda x: x.tag

    def __init__(self, path: str):
        self._tree    : etree._ElementTree = etree.parse(path, parser=None)
        self._head    : etree._Element     = self._tree.getroot()
        
        for child in self._head.getchildren():
            self.__dict__[child.tag] = XMLMetadataElement(child)

    def children(self):
        return list(map(XMLMetadata.__tag, self._head.getchildren()))

    def __getitem__(self, key: str):
        assert key in self.children(), f"Key <{key}> does not exist in tree."
        return XMLMetadataElement(self._head.xpath(key)[0])


class XMLMetadataElement(ABC):
    def __init__(self, element: etree._Element):
        self.element   = element
        self.text      = self.element.text
        self._children = list(map(lambda x: x.tag, self.element.getchildren()))
        __count        = [self._children.count(x)-1 for x in self._children]
        for i, child in enumerate(self.element.getchildren()):
            key = child.tag if not any(__count) else child.tag + '_%s' % i
            self.__dict__[key] = XMLMetadataElement(child)

    def __len__(self):
        return len(self._children)

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return f"<{type(self).__name__}>{self._children or self.element.text}"

    
class Input(ABC):
    "Base class for inputs."
    def __init__(self, PATH: str) -> None:
        pass
    
    @abstractmethod
    def __check(self):
        pass