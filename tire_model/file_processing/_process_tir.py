class _Processor:
    def __init__(self, name: str, file_path: str) -> None:
        self._tire = None
        self._exclude = ["$", "!"]

        self._add_tire(name, file_path)

    def _add_tire(self, name: str, path: str) -> None:
        self._tire = [name, self._import_data(path)]
    
    def get_parameters(self, parameter: str) -> dict:
        return self._tire[1][parameter]

    def _import_data(self, path: str) -> list[str, list]:
        local_results = {}
        f = open(path, "r")

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