from ROOT import * 
import pylab as plt

input_file = TFile('../output/output.root')  #reads in the root file with histograms
list_of_keys = input_file.GetListOfKeys()          #gets a list of the keys of all the objects in the file

#Creates an array containing the histogram objects and one of the radii. The Radius is taken from the last 4 string characters of the filename, so the Filenaming convention must be Histogram_name_blabla_x.xx where only x.xx is read out.
Histogramm_List = []
radii_list = []
for element in list_of_keys:
		Histogramm_List.append(input_file.Get(element.GetName()))
		radii_list.append(element.GetName()[-4:])

#Creates an array of the rms and mean data of all the objects
list_of_sigma_mean_touples = []
for element in Histogramm_List:
		list_of_sigma_mean_touples.append((element.GetStdDev(), element.GetMean()))

def sigma_over_mean(sigma_mean_touple):
		if sigma_mean_touple[1] == 0:
				result = 0
		else:
				result = sigma_mean_touple[0] / sigma_mean_touple[1]
		return result
#saves all the histograms, uncomment if necessary
'''
save = TFile("save.root", "new")	
for element in Histogramm_List:
		element.Write()
'''


#plots a graph of the uncertainty in the energy over the energy. 
sigma_over_mean_array = [sigma_over_mean(params) for params in list_of_sigma_mean_touples]
#comment out next 2 lines once files follow naming convention.
radii_listPlaceholder = range(len(radii_list))
plt.scatter(radii_listPlaceholder, sigma_over_mean_array)


#uncomment next 2 lines once files follow naming convention
#plt.scatter(radii_list, sigma_over_mean_array)
#radii_list = [float(radius_string) for radius_string in radii_list]

plt.ylim(min(sigma_over_mean_array), max(sigma_over_mean_array))
plt.xlim(0, max(radii_listPlaceholder))
plt.title("Sigma_E/E for various R")
plt.xlabel("radius (cm)")
plt.ylabel("Sigma_E/E")
plt.show()
