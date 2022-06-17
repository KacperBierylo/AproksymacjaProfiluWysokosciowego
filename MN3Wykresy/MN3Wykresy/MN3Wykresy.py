import matplotlib.pyplot as plt
import pandas as pd 
  
plikiOryginalne = ["SpacerniakGdansk", "kolorado2", "uluru"]
plikiWygenerowaneGdansk = [
   "GdanskLagrange5Wezlow", "GdanskLagrange5WezlowLos","GdanskLagrange10Wezlow","GdanskLagrange10WezlowLos",
   "GdanskLagrange15Wezlow", "GdanskLagrange15WezlowLos","GdanskLagrange20Wezlow", "GdanskLagrange20WezlowLos",

   "GdanskSklejanie5Wezlow", "GdanskSklejanie5WezlowLos","GdanskSklejanie10Wezlow","GdanskSklejanie10WezlowLos",
   "GdanskSklejanie15Wezlow", "GdanskSklejanie15WezlowLos","GdanskSklejanie20Wezlow", "GdanskSklejanie20WezlowLos",
]

plikiWygenerowaneKanion = [
   "KanionLagrange5Wezlow", "KanionLagrange5WezlowLos","KanionLagrange10Wezlow","KanionLagrange10WezlowLos",
   "KanionLagrange15Wezlow", "KanionLagrange15WezlowLos","KanionLagrange20Wezlow", "KanionLagrange20WezlowLos",

   "KanionSklejanie5Wezlow", "KanionSklejanie5WezlowLos","KanionSklejanie10Wezlow","KanionSklejanie10WezlowLos",
   "KanionSklejanie15Wezlow", "KanionSklejanie15WezlowLos","KanionSklejanie20Wezlow", "KanionSklejanie20WezlowLos",
]

plikiWygenerowaneUluru = [
   "UluruLagrange5Wezlow", "UluruLagrange5WezlowLos","UluruLagrange10Wezlow","UluruLagrange10WezlowLos",
   "UluruLagrange15Wezlow", "UluruLagrange15WezlowLos","UluruLagrange20Wezlow", "UluruLagrange20WezlowLos",

   "UluruSklejanie5Wezlow", "UluruSklejanie5WezlowLos","UluruSklejanie10Wezlow","UluruSklejanie10WezlowLos",
   "UluruSklejanie15Wezlow", "UluruSklejanie15WezlowLos","UluruSklejanie20Wezlow", "UluruSklejanie20WezlowLos",
]

for plik in plikiWygenerowaneGdansk:
    ramka = pd.read_csv("wyniki/" + plik + ".csv")  
    df = ramka
    x = df["Dystans (m)"].tolist()
    y = df["Wysokosc (m)"].tolist()
    xp = df["Dystans (m)"].tolist()
    yp = df["Wysokosc (m)"].tolist()
    xp2 = df["Dystans (m)"].tolist()
    yp2 = df["Wysokosc (m)"].tolist()

    wezly = df["Wezel"].tolist()

    j = 0;
    dlugosc = len(xp)
    for i in range(dlugosc):
        if wezly[j] == 0:
            xp.remove(xp2[j])
            yp.remove(yp2[j])
        j += 1

    plikorg = plikiOryginalne[0];
    ramka = pd.read_csv(plikorg + ".csv")  
    df = ramka
    x2 = df["Dystans (m)"].tolist()
    y2 = df["Wysokosc (m)"].tolist()
    fig2, ax = plt.subplots()
    ax.plot(x, y,"b", label="Wygenerowany")
    ax.plot(x2, y2,"r", label="Oryginalny")
    ax.plot(xp, yp, 'go')
    plt.xlabel('dystans [m]')
    plt.ylabel('wysokosc [m]')
    plt.title(plik)
    plt.savefig("wykresy/" + plik + '.png')


for plik in plikiWygenerowaneKanion:
    ramka = pd.read_csv("wyniki/" + plik + ".csv")  
    df = ramka
    x = df["Dystans (m)"].tolist()
    y = df["Wysokosc (m)"].tolist()
    xp = df["Dystans (m)"].tolist()
    yp = df["Wysokosc (m)"].tolist()
    xp2 = df["Dystans (m)"].tolist()
    yp2 = df["Wysokosc (m)"].tolist()

    wezly = df["Wezel"].tolist()

    j = 0;
    dlugosc = len(xp)
    for i in range(dlugosc):
        if wezly[j] == 0:
            xp.remove(xp2[j])
            yp.remove(yp2[j])
        j += 1

    plikorg = plikiOryginalne[1];
    ramka = pd.read_csv(plikorg + ".csv")  
    df = ramka
    x2 = df["Dystans (m)"].tolist()
    y2 = df["Wysokosc (m)"].tolist()
    fig2, ax = plt.subplots()
    ax.plot(x, y,"b", label="Wygenerowany")
    ax.plot(x2, y2,"r", label="Oryginalny")
    ax.plot(xp, yp, 'go')
    plt.xlabel('dystans [m]')
    plt.ylabel('wysokosc [m]')
    plt.title(plik)
    plt.savefig("wykresy/" + plik + '.png')

for plik in plikiWygenerowaneUluru:
    ramka = pd.read_csv("wyniki/" + plik + ".csv")  
    df = ramka
    x = df["Dystans (m)"].tolist()
    y = df["Wysokosc (m)"].tolist()
    xp = df["Dystans (m)"].tolist()
    yp = df["Wysokosc (m)"].tolist()
    xp2 = df["Dystans (m)"].tolist()
    yp2 = df["Wysokosc (m)"].tolist()

    wezly = df["Wezel"].tolist()

    j = 0;
    dlugosc = len(xp)
    for i in range(dlugosc):
        if wezly[j] == 0:
            xp.remove(xp2[j])
            yp.remove(yp2[j])
        j += 1

    plikorg = plikiOryginalne[2];
    ramka = pd.read_csv(plikorg + ".csv")  
    df = ramka
    x2 = df["Dystans (m)"].tolist()
    y2 = df["Wysokosc (m)"].tolist()
    fig2, ax = plt.subplots()
    ax.plot(x, y,"b", label="Wygenerowany")
    ax.plot(x2, y2,"r", label="Oryginalny")
    ax.plot(xp, yp, 'go')
    plt.xlabel('dystans [m]')
    plt.ylabel('wysokosc [m]')
    plt.title(plik)
    plt.savefig("wykresy/" + plik + '.png')
