from   abc  import ABC, abstractmethod
from   glob import glob
from   lxml import etree

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

    @property
    @abstractmethod
    def bands(self):
        ...


class XMLMetadata(ABC):
    "Base class for XML metadata. Base XML parser."
    __tag = lambda x: x.tag

    def __init__(self, path: str):
        self._tree    : etree._ElementTree = etree.parse(path)
        self._root    : etree._Element     = self._tree.getroot()
        self._children: etree._Element     = self.children()
        print(self._children)

    def children(self, key=None):
        element = self._root.find(key) if key else self._root
        return list(map(XMLMetadata.__tag, element.getchildren()))

    def __getitem__(self, *keys: str):
        assert keys in self.children, f"<{keys}> does not exist in tree."
        # return self._root.find('/'.join(key))
        return self._root.xpath('/'.join(key))





