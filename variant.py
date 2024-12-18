import numpy as np
import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt


# read_csv fonksiyonuyla birlikte "variants.tsv" dosyasini okudum.Okuma sonunda df degiskenine atadim.

df = pd.read_csv("variants.tsv", sep='\t')


                         ########## Veriye Hizli Bakis ############
# Bu fonksiyonu verinin genel resmine bakmak icin yazdim Cunku veriler uzerinde calisirken
# bazi degisikler yaptigimizda bu degisiklerinde acaba veri setinde bir degisiklige sebep olup olmadigini
#gozlemlemek istersem bu fonksiyonu calistirarask cok daha hizli bir sekilde verimi gozlemlemis olurum
def check_df(dataframe, head=5):
    print("######## Shape #############")
    print(dataframe.shape)
    print("########### Types ##############")
    print(dataframe.dtypes)
    print("############# Head ##############")
    print(dataframe.head(head))
    print("######### Tail ##################")
    print(dataframe.tail(head))
    print("############# NA #################")
    print(dataframe.isnull().sum())
    print("############# Quantiles ############")
    print(dataframe.describe([0, 0.05, 0.50, 0.95, 0.99, 1]).T)

check_df(df)

# GOREV1- Frequency sütunundaki eksik değerleri, o varyantın bulunduğu popülasyonın ortalama frekansı ile doldurun.


df[["Population", "Frequency"]].head() #gozlemlemek istedim

df[["Population", "Frequency"]].isnull().sum() #Eksik degerlerini gozlemlemek istedim

#gorupby fonksiyonuyla populasyonumu essiz bir sekilde gruplandirdim ardindan her grup için "Frequency" sütununu sectim
# transform fonksiyonuyla eksik frekans degerlerimi veri yapimi bozmadan dolduruyorum (apply kullansaydim veri yapim bozulurdu)
#transform(lambda x: x.fillna(x.mean())): Bu, her grup için eksik değerleri  grup ortalaması ile doldurur:
# lambda kullan at fonksiyonlarindan cok pratik fonksiyondur.Secilen populasyonun frekanslarinda gezin (x) gezdigin deger NaN deger ise onu ortalamayla doldur
# fillna fonksiyonu da NaN deger olup olmadigini kontrol ediyor

df["Frequency"] = df.groupby("Population")["Frequency"].transform(lambda x: x.fillna(x.mean()))


###### GOREV2- Clinical_Significance sütunundaki eksik değerleri, o popülasyondaki en sık görülen Clinical_Significance değeriyle doldurun.##########

#groupby fonksiyonuyla populasyonu gruplandirip "Clinical_Significance" sutununu sectim
# transform  fonkiyonunu her grup icin donusturmek adina kullandim
# lambda fonksiyonuyla her populasyonun sectigim clinical_signifance de geziniyorum ve sorguluyorum (fillna) eksik deger mi?
#eksik deger oldugunu anladiginda o sectigim x in yerine grubun mod degerini yerlestiriyor
#x.mode()[0] mod'un ilk öğesini seçer (çünkü mod, birden fazla öğe döndürebilir).

df["Clinical_Significance"].isnull().sum() # kac tane eksik degerinin oldugunu gozlemlemek istedim
df["Clinical_Significance"] = df.groupby("Population")["Clinical_Significance"].transform(lambda x: x.fillna(x.mode()[0]))
print(df)


####### GOREV3-Her popülasyonun ortalama frekans değerini hesaplayın. Ve bu ortalama değerin üzerinde ve eşit frekansa sahip olan varyantları filtreleyin.
#Yalnızca Pathogenic olarak işaretlenmiş varyantları içeren bir alt küme oluşturun.


#Her populasyonun ortalama frekansini aldim
P_F_mean = df.groupby("Population")["Frequency"].mean()
print(P_F_mean)

# Her satirda islem yapmak icin apply kullandim ve axis = 1 dedim ve lambdayla birlikte her bir ornekte gezindik (x)
# aldigim degeri P_F_mean degeriyle karsilastirdim sarti sagladigi takdirde filtered_variants degiskenine atadim

filtered_variants = df[df.apply(lambda x: x["Frequency"] >= P_F_mean[x["Population"]], axis=1)]

 ## baktim bu kac tane bulmus
print(filtered_variants.sum())
len(filtered_variants)

# Sadece "Pathogenic" olarak işaretlenmiş varyantları içeren bir alt küme oluşturma
pathogenic_variants = filtered_variants[filtered_variants["Clinical_Significance"] == "Pathogenic"]


print(pathogenic_variants)
len(pathogenic_variants)

#GOREV4 -Filtrelenmiş varyantlar üzerinde, bir Venn diyagramı ile varyantların popülasyonlar arasında ortaklık bilgisinin görselleştirilmesini sağlayın.


# Filtrelenmiş DataFrame örneği
filtered_variants = pathogenic_variants  # Daha önce oluşturduğumuz patojenik varyantların alt kümesi

# Her popülasyondaki varyant setini oluşturdum
european_variants = set(filtered_variants[filtered_variants["Population"] == "European"]["Variant"])
african_variants = set(filtered_variants[filtered_variants["Population"] == "African"]["Variant"])
asian_variants = set(filtered_variants[filtered_variants["Population"] == "Asian"]["Variant"])

# Venn diyagramı oluşturdum
plt.figure(figsize=(8, 8))
venn = venn3(
    subsets=(european_variants, african_variants, asian_variants),
    set_labels=("European", "African", "Asian"))


# Venn diagramını çizmek için matplotlib_venn modülünü kullandım.
plt.figure(figsize=(10, 8))

# Üç popülasyonun varyant setlerini bir Venn diagramında gösteriyordum.
venn = venn3([european_variants, african_variants, asian_variants],
             set_labels=("European", "African", "Asian"))

# Başlık ekledım
plt.title("Varyantların Popülasyonlar Arasındaki Ortaklığı")


plt.show()

##########################################



