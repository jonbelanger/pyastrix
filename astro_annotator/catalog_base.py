from abc import ABC, abstractmethod


class CatalogQuery(ABC):
    """
    Standardized output fields REQUIRED:
        - source_id
        - ra (deg)
        - dec (deg)
    """

    @abstractmethod
    def query(self, center, radius):
        """
        Must return a table-like object with
        standardized columns: source_id, ra, dec
        """
        pass
