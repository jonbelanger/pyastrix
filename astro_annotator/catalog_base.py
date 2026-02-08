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

    def format_label(self, row):
        """Optional: produce a single display label for a result row.

        Catalog implementations can override this to return a string.
        If `None` is returned, callers may fall back to other formatting.
        """
        return None
