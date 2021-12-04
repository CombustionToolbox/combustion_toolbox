function index = get_index_ions(self)
    index = (contains(self.S.LS, 'minus') | contains(self.S.LS, 'plus')) & ~contains(self.S.LS, 'cyclominus');
end