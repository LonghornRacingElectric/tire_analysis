from typing import Sequence, Tuple


class Processor:
    """
    ## Processor

    Processes .tir files

    Parameters
    ----------
    name : str
        Internal name of processed result (.tir values)
    file_path : str
        Path to desired .tir file
    """
    def __init__(self, name: str, file_path: str) -> None:
        self._tire = None
        self._exclude = ["$", "!"]
        self.file_path = file_path

        self._add_tire(name)

    def _add_tire(self, name: str) -> None:
        """
        ## Add Tire

        Stores desired tire locally

        Parameters
        ----------
        name : str
            Internal name of processed result
        """
        self._tire = [name, self._import_data()]
    
    def get_parameters(self, parameter: str) -> dict:
        """
        ## Get Parameter

        Displays desired parameter from tire properties

        Parameters
        ----------
        parameter : str
            Parameter to display

        Returns
        -------
        dict
            Dictionary of (key, value) pairs for desired parameter
        """
        return self._tire[1][parameter]

    def _import_data(self) -> Sequence[Tuple[str, list]]:
        """
        ## Import Data

        Parses .tir file and stores parameters

        Parameters
        ----------
        None

        Returns
        -------
        Sequence[Sequence[str, list]]
            Sequence containing .tir headers and corresponding (parameter, value) pairs
        """
        local_results = {}
        f = open(self.file_path, "r")

        data_entry = False

        for line in f:
            char_0 = line.strip()[0]

            if data_entry and (char_0 not in self._exclude):
                line_stripped = line.replace(" ", "")

                if "$" in line_stripped:
                    line_stripped = line_stripped[:line_stripped.index("$")]

                line_split = line_stripped.split("=")
                if line_split[1].replace(".", "").replace("-", "").replace("E", "").replace("e", "").replace("+", "").replace("\n", "").isnumeric():
                    val = float(line_split[1])
                
                else:
                    val = line_split[1]

                local_results[list(local_results.keys())[-1]][line_split[0]] = val

            else:
                if (char_0 in self._exclude):
                    data_entry = False
                    continue
            
                if (char_0 == "["):
                    if ("[SHAPE]" in line):
                        continue

                    local_results[line.strip()[1:-1]] = {}
                    data_entry = True
                
        return local_results