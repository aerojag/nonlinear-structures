#

def trial_x(self, i, j):
    return (1 - self.x) ** i * (1 + self.x) ** j * self.P[self.i]
    # if i == 0 and j == 0:
    #     return self.P[self.i]
    # elif i == 0:
    #     return (1 + self.x) ** j * self.P[self.i]
    # elif j == 0:
    #     return (1 - self.x) ** i * self.P[self.i]
    # else:
    #     return (1 - self.x) ** i * (1 + self.x) ** j * self.P[self.i]

def trial_dx(self, i, j):
    return - i * (1 - self.x) ** (i - 1) * (1 + self.x) ** j * self.P[self.i] + \
                j * (1 - self.x) ** i * (1 + self.x) ** (j - 1) * self.P[self.i] + \
                (1 - self.x) ** i * (1 + self.x) ** j * self.dP[self.i]

def trial_dx(self, i, j):
        return - i * (1 - self.x) ** (i - 1) * (1 + self.x) ** j * self.P[self.i] if i else 0 + \
                 j * (1 - self.x) ** i * (1 + self.x) ** (j - 1) * self.P[self.i] if j else 0 + \
                ((1 - self.x) ** i if i else 1) * ((1 + self.x) ** j if j else 1) * self.dP[self.i]
        # if i == 0 and j == 0:
        #     return self.dP[self.i]
        # elif i == 0:
        #     return self.P[self.i] + (1 + self.x) * self.dP[self.i]
        # elif j == 0:
        #     return - self.P[self.i] + (1 - self.x) * self.dP[self.i]
        # else:
        #     return - 2 * self.x * self.P[self.i] + (1 - self.x) * (1 + self.x) * self.dP[self.i]


def trial_d2x(self, i, j):
    return i * (i - 1) * (1 - self.x) * (i - 2) * (1 + x) * j * self.P[self.i] + \
            j * (j - 1) * (1 - self.x) * i * (1 + self.x) * (j - 2) * self.P[self.i] - \
            2 * i * (1 - self.x) * (i - 1) * (1 + self.x) * j * self.dP[self.i] + \
            2 * j * (1 - self.x) * i * (1 + self.x) * (j - 1) * self.dP[self.i] - \
            2 * i * j * (1 - self.x) * (i - 1) * (1 + self.x) * (j - 1) * self.P[self.i] + \
            (1 - self.x) * i * (1 + self.x) * j * self.d2P[self.i]
    