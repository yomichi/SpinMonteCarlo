class Square:
    def __init__(self, L):
        self.L = L
        self.N = L*L

    def num_site(self) : return self.N
    def num_bond(self) : return 2*self.N

    def coord2site(self, x,y):
        return x + self.L*y

    def site2coord(self, site):
        return site%self.L, site/self.L

    def neighbors(self, site):
        x,y = self.site2coord(site)
        return  [ self.coord2site( (x+self.L-1)%self.L, y),
                  self.coord2site( (x+1)%self.L, y),
                  self.coord2site( x, (y+self.L-1)%self.L),
                  self.coord2site( x, (y+1)%self.L)]
