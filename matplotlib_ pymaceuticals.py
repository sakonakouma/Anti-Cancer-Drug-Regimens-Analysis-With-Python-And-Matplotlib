#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import dependencies 
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import os


# In[2]:


# OBSERVATIONS
#Look across all previously generated figures and tables and write at least three observations or inferences 
#that can be made from the data. Include these observations at the top of notebook.

# Observation 1: The data analysis shows positive relationship between tumor volume and weights. Indeed, the the 
# coefficient of correlation between the mouse weight and the average tumor volume  (R2) = 0.84

# Observation 2: Pvalue = 1.32 which shows the model do not explain the relationship between the mouse weight and 
# the average tumor volume.

# Observation 3: The plot line shows their is no relationship between the timepoint and the tumor volume. 
# Also the box plot and whisker shows that only one drug regimen (infubinol) has outliers or extreme values.
# Also, the data analysis shows the mouse distribution by sex is almost equal.


# In[3]:


os.getcwd()


# In[4]:


# Set the path to the file
os.chdir("C:\\Users\\sakon\\Desktop\\MASTER CSU\\CWRU_BOOTCAMP\\1.Hmw\Matplotlib-challenge\\")


# In[5]:


file1=os.path.join("C:\\Users\\sakon\\Desktop\\MASTER CSU\\CWRU_BOOTCAMP\\1.Hmw\Matplotlib-challenge","Study_results.csv")


# In[6]:


file2=os.path.join("C:\\Users\\sakon\\Desktop\\MASTER CSU\\CWRU_BOOTCAMP\\1.Hmw\Matplotlib-challenge","Mouse_metadata.csv")


# In[7]:


print (file1)


# In[8]:


Study_results = pd.read_csv("C:\\Users\\sakon\\Desktop\\MASTER CSU\\CWRU_BOOTCAMP\\1.Hmw\Matplotlib-challenge\\Study_results.csv")
Mouse_metadata = pd.read_csv("C:\\Users\\sakon\\Desktop\\MASTER CSU\\CWRU_BOOTCAMP\\1.Hmw\Matplotlib-challenge\\Mouse_metadata.csv")


# In[9]:


Study_results.head()


# In[10]:


Mouse_metadata.head()


# In[11]:


#Merge Mouse metadata and study results

DataMerge1=pd.merge(Mouse_metadata, Study_results, how='right', on=['Mouse ID'])
DataMerge1


# In[12]:


Data=DataMerge1.drop_duplicates()


# In[13]:


Data.shape


# In[14]:


# Group data by drug regimen

Datagroup=Data.groupby("Drug Regimen")
Datagroup


# In[15]:


DataStat1=Datagroup["Tumor Volume (mm3)"].sum()
DataStat2=Datagroup["Tumor Volume (mm3)"].mean()
DataStat3=Datagroup["Tumor Volume (mm3)"].median()
DataStat4=Datagroup["Tumor Volume (mm3)"].var()
DataStat5=Datagroup["Tumor Volume (mm3)"].sem()
DataStat6=Datagroup["Tumor Volume (mm3)"].std()


# In[16]:


print(DataStat1)


# In[17]:


#Generate a summary statistics table consisting of the mean, median, variance, standard deviation, 
# and SEM of the tumor volume for each drug regimen

Summary_table=pd.DataFrame({"Total tumor volume":DataStat1, 
                            "Tumor Mean":DataStat2, 
                            "Tumor Median":DataStat3, 
                            "Tumor Variance":DataStat4,
                            "Tumor standard deviation":DataStat6,
                            "Tumor Standard Error of the Mean": DataStat5})
round(Summary_table, 2)


# In[18]:


DataStat7=Datagroup["Tumor Volume (mm3)"].count()
DataStat7


# In[19]:


Myplot=DataStat7.plot.bar(color='Red')
Myplot
plt.ylabel("Number of data points")
plt.title("Data Points per Drug Regiment")
plt.savefig('barplot1')


# In[20]:


DataStat7.plot.bar(color='Blue', title='Data Points per Drug Regiment', legend=True, grid=True,)


# In[21]:


Data_pie=Mouse_metadata.groupby('Sex').count()
Data_pie


# In[22]:


labels=[Mouse_metadata['Sex'].unique()]
L1=list(Data_pie.index)


# In[23]:


Sizes=[Data_pie['Mouse ID']]
Sizes


# In[24]:


Colors=['orange','blue']
plt.pie(Sizes, labels=L1, colors=Colors, autopct='%1.1f%%', shadow=True, startangle=180)
plt.title('Distribution by sex in the mouse population')
plt.ylabel('Sex')
plt.show


# In[25]:


Matplot_pie=Data_pie['Mouse ID'].plot(kind='pie', y='Sizes', figsize=(6,6), autopct='%1.1f%%', shadow=True, startangle=180)
plt.title('Distribution by sex in the mouse population')
plt.ylabel('Sex')
plt.show


# In[26]:


DataPro=Data[Data['Drug Regimen'].isin(['Capomulin', 'Ramicane', 'Infubinol',
                                        'Ceftamin'])].sort_values(['Timepoint'],ascending=True)
DataPro


# In[27]:


MouseWei=DataPro.loc[DataPro['Drug Regimen']=="Capomulin",:]
MouseWei 


# In[28]:


# Calculate the quartiles and IQR and quantitatively determine if there are any potential outliers across all four treatment 
# regimens.
MouseWeigrpe=MouseWei.groupby("Mouse ID").max()["Timepoint"]
MouseWeigrpe_df=pd.DataFrame(MouseWeigrpe)
MouseWeigrpe_df


# In[29]:


MouseScatter=pd.merge(MouseWeigrpe_df, MouseWei, on=('Mouse ID', 'Timepoint'))
MouseScatter


# In[30]:


#Calculate the quartiles and IQR and quantitatively determine if there are any potential outliers across all 
# four treatment regimens.
StatQuart= MouseScatter['Tumor Volume (mm3)']
Quartiles= StatQuart.quantile([0.25, 0.50, 0.75])
LowerQ=Quartiles[0.25]
UpperQ=Quartiles[0.75]
IQR1=UpperQ-LowerQ
Lowerbound=LowerQ-(1.5*IQR1)
Upperbound=UpperQ+(1.5*IQR1)
print(f"Capomulin outliers values below{Lowerbound} and above {Upperbound}")


# In[31]:


MouseWeiR=DataPro.loc[DataPro['Drug Regimen']=="Ramicane",:]
MouseWeiR


# In[32]:


MouseWeigrpeR=MouseWeiR.groupby("Mouse ID").max()["Timepoint"]
MouseWeigrpe_dfR=pd.DataFrame(MouseWeigrpeR)
MouseWeigrpe_dfR


# In[33]:


MouseScatterR=pd.merge(MouseWeigrpe_dfR, MouseWeiR, on=('Mouse ID', 'Timepoint'))
MouseScatterR


# In[34]:


StatQuartR= MouseScatterR['Tumor Volume (mm3)']
QuartilesR= StatQuartR.quantile([0.25, 0.50, 0.75])
LowerQR=QuartilesR[0.25]
UpperQR=QuartilesR[0.75]
IQR1R=UpperQR-LowerQR
Lowerbound=LowerQR-(1.5*IQR1R)
Upperbound=UpperQR+(1.5*IQR1R)
print(f"Ramicane outliers values below {Lowerbound} and above {Upperbound}")


# In[35]:


MouseWeiI=DataPro.loc[DataPro['Drug Regimen']=="Infubinol",:]
MouseWeiI 


# In[36]:


MouseWeigrpeI=MouseWeiI.groupby("Mouse ID").max()["Timepoint"]
MouseWeigrpe_dfI=pd.DataFrame(MouseWeigrpeI)
MouseWeigrpe_dfI


# In[37]:


MouseScatterI=pd.merge(MouseWeigrpe_dfI, MouseWeiI, on=('Mouse ID', 'Timepoint'))
MouseScatterI


# In[38]:


StatQuartI= MouseScatterI['Tumor Volume (mm3)']
Quartiles= StatQuartI.quantile([0.25, 0.50, 0.75])
LowerQ=Quartiles[0.25]
UpperQ=Quartiles[0.75]
IQR1=UpperQ-LowerQ
Lowerbound=LowerQ-(1.5*IQR1)
Upperbound=UpperQ+(1.5*IQR1)
print(f"Infubinol outliers values below{Lowerbound} and above {Upperbound}")


# In[39]:


MouseWeiC=DataPro.loc[DataPro['Drug Regimen']=="Ceftamin",:]
MouseWeiC 


# In[40]:


MouseWeigrpeC=MouseWeiC.groupby("Mouse ID").max()["Timepoint"]
MouseWeigrpe_dfC=pd.DataFrame(MouseWeigrpeC)
MouseWeigrpe_dfC


# In[41]:


MouseScatterC=pd.merge(MouseWeigrpe_dfC, MouseWeiC, on=('Mouse ID', 'Timepoint'))
MouseScatterC


# In[42]:


StatQuartC= MouseScatterC['Tumor Volume (mm3)']
Quartiles= StatQuartC.quantile([0.25, 0.50, 0.75])
LowerQ=Quartiles[0.25]
UpperQ=Quartiles[0.75]
IQR1=UpperQ-LowerQ
Lowerbound=LowerQ-(1.5*IQR1)
Upperbound=UpperQ+(1.5*IQR1)
print(f"Ceftamin outliers values below {Lowerbound} and above {Upperbound}")


# In[43]:


DataPro1=DataPro[['Mouse ID', 'Drug Regimen', 'Timepoint', 'Tumor Volume (mm3)']]
DataPro1


# In[ ]:





# In[44]:


Data1_sort=DataPro1.groupby(["Drug Regimen", "Mouse ID"]).last()['Tumor Volume (mm3)']
Data1_sort


# In[45]:


Data1_sort_df=Data1_sort.to_frame()
Data1_sort_df


# In[46]:


# generate a box and whisker plot of the final tumor volume for all four treatment regimens and highlight 
# any potential outliers 
name=['Capomulin','Ramicane', 'Infubinol', 'Ceftamin']
DataF=Data1_sort_df.reset_index()
Data_list=Data1_sort_df.groupby('Drug Regimen')['Tumor Volume (mm3)'].apply(list)
Data_list_df=pd.DataFrame(Data_list)
Data_list_df=Data_list_df.reindex(name)
DataFinal=[x for x in Data_list_df['Tumor Volume (mm3)']]
plt.boxplot(DataFinal, labels=name)
plt.ylim(10,80)
plt.show()


# In[47]:


#Generate a line plot of time point versus tumor volume for a single mouse treated with Capomulin. 

M957=DataPro1.loc[DataPro1['Mouse ID']=="m957",:]
M957


# In[48]:


#Generate a line plot of time point versus tumor volume for a single mouse treated with Capomulin.

M957.plot.line(x='Timepoint',y='Tumor Volume (mm3)')


# In[49]:


# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin treatment regimen.

MouseWei=DataPro.loc[DataPro['Drug Regimen']=="Capomulin",:]
MouseWei


# In[50]:


#Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin treatment regimen.
DataScart=MouseWei.groupby(["Mouse ID"]).mean()
plt.scatter(DataScart['Weight (g)'], DataScart['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Tumor Volume (mm3)')


# In[51]:


# Calculate the correlation coefficient and linear regression model between mouse weight and average tumor volume for the 
# Capomulin treatment. Plot the linear regression model on top of the previous scatter plot.
import scipy.stats as st
Cor=st.pearsonr(DataScart['Weight (g)'], DataScart['Tumor Volume (mm3)'])[0]
print(f"the coefficient of correlation between the mouse weight and the average tumor volume is {Cor}")


# In[52]:


LinReg=st.linregress((DataScart['Weight (g)'], DataScart['Tumor Volume (mm3)']))
LinReg


# In[53]:


#Calculate the correlation coefficient and linear regression model between mouse weight and average tumor volume for the Capomulin treatment. Plot the linear regression model on top of the previous 
# scatter plot.a=0.9544396890241045
b=21.552160532685015
y=DataScart['Weight (g)']*a+b
plt.scatter(DataScart['Weight (g)'], DataScart['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')
plt.plot(DataScart['Weight (g)'], y, color='red')
plt.show()

