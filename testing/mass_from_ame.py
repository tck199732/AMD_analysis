

from pyamd.utilities import ame


def main():
    mass_table = ame.AME()

    mass_4He = mass_table.get_mass('He4')
    mass_6Li = mass_table.get_mass('Li6')
    mass_p = mass_table.get_mass('H1')
    mass_d = mass_table.get_mass('H2')
    mass_t = mass_table.get_mass('H3')
    mass_3He = mass_table.get_mass('He3')
    mass_n = mass_table.get_mass('n1')

    mass_8Be = mass_table.get_mass('Be8')
    mass_7Be = mass_table.get_mass('Be7')
    mass_7Li = mass_table.get_mass('Li7')
    mass_6Li = mass_table.get_mass('Li6')
    mass_5Li = mass_table.get_mass('Li5')
    mass_12C = mass_table.get_mass('C12')
    print(mass_n)
    print(mass_p)
    print(mass_d)
    print(mass_t)
    print(mass_3He)
    print(mass_4He)

    # print(mass_5Li)
    # print(mass_6Li)
    # print(mass_7Li)
    # print(mass_7Be)
    # print(mass_8Be)
    print(mass_12C)


if __name__ == '__main__':
    main()
