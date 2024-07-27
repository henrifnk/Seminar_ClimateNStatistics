---
title: "Climate And Statistics"
author: "Helmut Küchenhoff, Henri Funk"
date: "2024-07-19"
documentclass: krantz
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
colorlinks: yes
lot: False
lof: False
site: bookdown::bookdown_site
description: "A Seminar about statistical methods in climate research in SS24."
graphics: yes
---
<!--- cover-image: images/cover.png -->



# Preface {-}

*Author: Henri Funk*


\begin{center}\includegraphics[width=500]{cover} \end{center}

As the world faces the reality of climate change, natural hazards and extreme weather events have become a major concern, with devastating consequences for nature and humans. The quantification and definition of climate change, extreme events and its implications for life and health on our planet is one of the major concerns in climate science. 

This book explains current statistical methods in climate science and their application.
The methods include compound events, low flow events and return periods, natural variability, teleconnections and causal discovery.
All of those methods are used to quantify and anticipate the changing climate.

This book is the outcome of the seminar "Climate and Statistics" which took place in summer 2024 at the Department of Statistics, LMU Munich.

![Creative Commons License](by-nc-sa.png)

This book is licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).


\mainmatter

# Foreword {-}

*Author: Christoph Molnar*

<!-- An experiment -->
This book is the result of an experiment in university teaching.
Each semester, students of the Statistics Master can choose from a selection of seminar topics.
Usually, every student in the seminar chooses a scientific paper, gives a talk about the paper and summarizes it in the form of a seminar paper.
The supervisors help the students, they listen to the talks, read the seminar papers, grade the work and then ... hide the seminar papers away in (digital) drawers.
This seemed wasteful to us, given the huge amount of effort the students usually invest in seminars.
An idea was born:
Why not create a book with a website as the outcome of the seminar?
Something that will last at least a few years after the end of the semester.
In the summer term 2019, some Statistics Master students signed up for our seminar entitled "Limitations of Interpretable Machine Learning".
When they came to the kick-off meeting, they had no idea that they would write a book by the end of the semester.

We were bound by the examination rules for conducting the seminar, but otherwise we could deviate from the traditional format.
We deviated in several ways:

1. Each student project is part of a book, and not an isolated seminar paper.
1. We gave challenges to the students, instead of papers. The challenge was to investigate a specific limitation of interpretable machine learning methods.
1. We designed the work to live beyond the seminar.
1. We emphasized collaboration. Students wrote some chapters in teams and reviewed each others texts.

<!-- Our experience -->
<!---
Looking back, the seminar was a lot of fun and -- from our perspective -- successful.
Especially considering that it was an experiment.
Everyone was highly motivated and we got great feedback from the students that they liked the format.
For the students it was a more work than a traditional seminar.
But in the end, our hope is that their effort will pay off for them as well, not only because of their increased visibility.
It was also more work for us supervisors.
But the extra effort was worth it, since limitations of interpretability are relevant for our research.
For me the seminar was an inspiration.
The students had new ideas and new perspectives to approach the limitations of interpretable machine learning.
-->

<!-- Technical setup -->
## Technical Setup {-}

The book chapters are written in the Markdown language.
The simulations, data examples and visualizations were created with R [@rlang].
To combine R-code and Markdown, we used rmarkdown.
The book was compiled with the bookdown package.
We collaborated using git and github.
For details, head over to the [book's repository](link/to/repo).



<!--chapter:end:index.Rmd-->

# Introduction

*Author: *

*Supervisor: *

## Intro About the Seminar Topic


## Outline of the Booklet


<!--chapter:end:00-introduction.Rmd-->

# Natural Variability by internal variability {#iv}

*Author: Author*

*Supervisor: Henri Funk*

*Suggested degree: Bachelor*

Natural variability refers to the inherent fluctuations in the climate system that occur without external forcings, such as changes in solar radiation, volcanic eruptions, or human-induced alterations of the Earth's atmosphere and surface. This variability can be due to a variety of factors, including atmospheric processes, ocean currents, the El Niño-Southern Oscillation (ENSO), and other dynamic components of the Earth system. Natural variability occurs across all time scales, from short-term (daily, seasonal) to long-term (decadal, centennial) fluctuations.

## Climate Model Ensemble

Climate models are sophisticated tools that simulate the interactions within the Earth's climate system. To understand and quantify natural variability, scientists use ensembles of climate model simulations. An ensemble consists of multiple runs of the same model, or different models, where each run has slightly different initial conditions or model parameters. This approach helps to capture the range of possible climate outcomes due to the inherent uncertainty and variability in the system.

## Internal Variability
Within the context of climate model ensembles, internal variability refers to the variations in climate that arise from the system's internal processes, independent of external forcings. This variability is a fundamental part of the climate system's dynamics and can lead to different outcomes even if the external conditions are the same.

@deser


<!--chapter:end:01-variability.Rmd-->

# Standard Precipitation Evapotranspiration Index {#spei}

*Author: Author*

*Supervisor: Henri Funk*

*Suggested degree: Bachelor*

## Definition

The SPEI is a statistical indicator that shows dry and wet periods on the basis of the climatic water balance (total precipitation minus total potential evaporation). It captures water deficits or surpluses on the land surface better than the Standardised Precipitation Index (SPI), which focuses solely on precipitation. The potential evaporation can be based on the so-called FAO grass evaporation, which takes into account radiation, air temperature, relative humidity and wind speed over a standard grass area. The time series is broken down into 12 monthly frequency sums, to each of which a cumulative log-logistic distribution is fitted. This is transformed into a corresponding cumulative standard normal distribution, the abscissa value of which - referred to as SPEI - allows simple assignment to probability classes.
@vicente

<!--chapter:end:02-spei.Rmd-->

# Compound events {#ce}

*Author: Author*

*Supervisor: Henri Funk*

*Suggested degree: Master*

In climate research, a compound event refers to the combination of multiple extreme weather or climate events occurring simultaneously or successively, leading to significant impacts on the environment, society, or economy. These events can be related either through meteorological factors (e.g., heatwaves and drought occurring together due to a prolonged high-pressure system) or through their impacts (e.g., heavy rainfall and storm surge combining to cause more severe flooding). The complexity of compound events lies in their interconnectedness and the way they can exacerbate each other's effects, often in a nonlinear manner, making them more challenging to predict and manage than individual extreme events.

## Example: Hot and Dry Events

Consider a scenario where a region experiences both extreme heat and an extended period of low precipitation simultaneously. In this example, the compound nature of the hot and dry events creates a cascade of impacts that are more severe than what would be expected if these events occurred independently. The interplay between the heat and drought amplifies the consequences, underscoring the importance of understanding and managing compound events in the context of climate change and variability.
@zscheischler


<!--chapter:end:03-compounds.Rmd-->

# Assymetric bivariate Copulas for compounds {#ac}

*Author: Author*

*Supervisor: Henri Funk*

*Suggested degree: Master*

In compounds, particularly in environmental and hydrological contexts, bivariate relationships can be used to study how two different environmental factors, such as temperature and precipitation, interact with each other. These relationships are critical for modeling the joint behavior of these variables, which can be essential for predicting weather events, designing structures, or managing natural resources.

In the study of environmental and hydrological systems, recognizing and accurately modeling the asymmetry in bivariate relationships is crucial for making reliable predictions and informed decisions. Techniques such as copulas, are often used in this context to model and analyze these complex and dependencies effectively. To model the asymmetry in dependencies, Archimax could be used.

@charpentier
@bacigal

<!--chapter:end:04-archimax.Rmd-->

# Teleconnections North Atlantic Oscillation {#nao}

*Author: Author*

*Supervisor: Henri Funk*

*Suggested degree: Bachelor/Master*

The concept of teleconnections in climate science refers to climate anomalies or patterns that are related across large distances, often thousands of kilometers apart. Teleconnections represent atmospheric interactions that link weather and climate conditions in one region of the globe with those in another, often through the movement and behavior of large-scale atmospheric wave patterns. These connections are crucial for understanding regional climate variations, predicting weather patterns, and assessing climate change impacts.

## North Atlantic Oscillation (NAO) and Its Implications

The North Atlantic Oscillation (NAO) is a prime example of a teleconnection pattern, characterized by variations in the difference of atmospheric pressure at sea level between the Icelandic Low and the Azores High. The NAO influences weather and climate conditions across the North Atlantic region and beyond, affecting winter temperatures and storm tracks.

## Implications for Our Climate

Understanding the NAO and other teleconnection patterns is essential for accurate weather forecasting, climate prediction, and the development of strategies to mitigate the impacts of climate variability and change. Researchers use observations, climate models, and paleoclimate reconstructions to study these patterns, aiming to improve our ability to predict their future behavior under different climate change scenarios.
@hurrell

<!--chapter:end:05-nao.Rmd-->

# Low Flow Events {#lfe}

*Author: Author*

*Supervisor: Henri Funk*

*Suggested degree: Bachelor*

Low flow events in hydrology refer to periods when water flow in rivers, streams, or creeks falls below a critical threshold, often leading to various ecological, economic, and social impacts. These events can result from prolonged periods of below-average precipitation, increased water usage, or changes in land use that affect the hydrological cycle. Low flow conditions can compromise water quality, reduce water availability for agricultural, industrial, and domestic uses, and disrupt aquatic ecosystems.

## Concept of Extremeness by Return Periods

The "extremeness" of hydrological events, including low flow conditions, is often characterized by return periods. A return period, also known as a recurrence interval, is a statistical measure used to estimate the frequency at which a certain intensity of an event (e.g., flow rate, rainfall amount) is likely to be equaled or exceeded. It is typically expressed in years.

- **Definition**: The return period is calculated based on historical data and is intended to provide an estimate of the likelihood of experiencing an event of a certain magnitude or greater within any given year. For example, a 100-year return period for a low flow event means that, on average, the flow rate is expected to fall to that level or lower once every 100 years. This does not imply the events will occur exactly 100 years apart but rather conveys a 1% chance of occurrence in any given year.
@du

<!--chapter:end:06-lowflow.Rmd-->

# Statistical streamflow modelling {#sm}

*Author: Author*

*Supervisor: Henri Funk*

*Suggested degree: Master*

Hydrological models are tools used to understand, simulate, and predict the movement, distribution, and quality of water within natural and engineered water systems. These models help in water resource management, flood forecasting, environmental protection, and understanding the impact of climate change on water resources.

There are mainly two types of hydrological models: physical and statistical.

## Physical Hydrological Models

Physical hydrological models are based on the physical processes that occur within the hydrological cycle. They use mathematical representations of these processes to simulate the behavior of water in the environment. Physical models consider processes such as precipitation, evaporation, infiltration, surface runoff, and groundwater flow. These models are often more complex and require detailed data about the terrain, soil properties, vegetation, and meteorological conditions.

## Statistical Hydrological Models

Statistical hydrological models use historical data to identify patterns and relationships between different hydrological variables. These models apply statistical methods to predict future water conditions based on past observations. They do not necessarily simulate the physical processes but rather rely on the statistical properties of the data, such as trends, cycles, and variability. Statistical models are often used when the physical processes are too complex to model directly or when there is insufficient data to build a detailed physical model.

@sabzipour

<!--chapter:end:07-hydroLSTM.Rmd-->

# The Lancet Report 2023    {#he1}

*Author: Author*

*Supervisor: Helmut Kuechenhoff*

*Suggested degree: Bachelor* 

In yearly reports one of the leading scientific journals in health research 
gives an overview on the state of the art knowledge on the effect of climate 
change on health (@romanello).  The following aspects should be treated: 

#. Main results and message 
#.  Communication of results 
#. Arguments based on statistical methods 


<!--chapter:end:08-lancet.Rmd-->

#  Epidemiologic studies on the heat effects     {#he2}

*Author: Author*

*Supervisor: Helmut Kuechenhoff*

*Degree: Bachelor*   

There are many epidemiological studies for assessing the effect of heat on health.
One strategy is to use daily mortality data and weather data in certain areas to 
estimate heat related mortality by regression models. Two recent papers by groups of 
researchers give an overview: 

#. General effects 
@masselot

#. Joint effects of air pollution and heat 
@stafoggia
  
Results and statistical methods should be presented

<!--chapter:end:09-epidemiologic.Rmd-->

---
output:
  pdf_document: default
  word_document: default
  html_document: default
---
# Controversial issue : heat and humidity    {#he3}

*Author: Mona Niethammer*

*Supervisor: Helmut Kuechenhoff*

*Degree: Master*  



## Abstract 

As climate change progresses, heat-related mortality and adverse health outcomes are also increasing. In the chapters before it was mainly talked about temperature related outcomes. This chapter focuses on heat events and how humidity plays a crucial role in high ambient temperatures and the effect on health outcomes. Physiological research indicates that both heat and humidity significantly impact the body’s response to heat. Thus, influencing mortality and adverse health outcomes. Several epidemiological studies aimed to show this relationship, found however only minimal or even no effect of humidity on negative health outcomes. This discrepancy between physiological knowledge and epidemiological findings creates a, as (@baldwin2023) describes it, controversial issue. This will be further explored in this chapter.

## Introduction 

### Background 

Climate change is a broad and pressing topic as the number of heat-related deaths and the air-temperature break records. Extreme heat waves have become a consistent feature of summer seasons, leading to an increase in heat-related fatalities (@ebi2021). Rising greenhouse gas concentrations further contribute to rising temperatures and humidity at the same time. Physiological knowledge indicates that heat and humidity both contribute to human heat stress, adversely affecting the human body (@baldwin2023). Heat stress occurs when environmental conditions overwhelm the body‘s cooling mechanism (@buzan2020). When the body is not able to sufficiently lower its core temperature, the human body suffers health problems, potentially leading to death. From a physiological perspective, the body’s cooling mechanisms are clearly influenced both by heat and humidity. 

In contrast, a wide range of epidemiological studies have found either no or only a weak correlation between humidity and human heat stress or adverse health outcomes. Those studies all infer that a rise in air-temperature (heat) is associated with negative health outcomes. These two contradictions lead, according to (@baldwin2023), to a controversial issue as physiologically it is known that there is an effect of humidity on health outcomes, whereas in epidemiological studies they did not find evidence for such an effect. Figure 1 describes this contradiction between physiological knowledge and epidemiological studies. Heat in this chapter is used as synonym for high ambient temperature. 


\begin{center}\includegraphics[width=500]{Controversial Issue} \end{center}

**Figure 1:** <font size="2">Flowchart illustrating the controversial issue between epidemiological studies and physiological knowledge. \
*A: illustrates how heat and humidity changes with climate change. B: illustrates that epidemiological studies only take heat as driver into account, the physiological view heat and humidity. C: shows the predictions from an epidemiological and physiological view.* </font> 

### Epidemiological Studies 
Numerous epidemiological studies investigated the influence of heat and humidity on negative health outcomes such as death or cardiovascular diseases. Results of these studies are mixed regarding the association of humidity and negative health outcomes. As an example, a study in New York, conducted by Barreca found a positive correlation between high levels of relative humidity (RH) and cardio-respiratory hospitalizations. However, this study is one of few that identified an effect of humidity on negative health outcomes. 
Baldwin et al. reviewed several epidemiological studies and conducted further research. They concluded six key reasons why epidemiological studies found less or even no effect of humidity on health outcomes:

1) At high temperatures, humidity has minimal or even no influence on negative health outcomes \
2) Study results may be skewed by limited data sets \
3) Analyses often focus on vulnerable populations, like elderly, infants or individuals with impaired sweating functions \
4) Necessary extreme values of heat and humidity to show a significant impact on heat strain are rare \
5) The relationship between heat and humidity is often not adequately considered. This may lead to inappropriate model results \
6) Sub-daily meteorological phenomena, such as rain, which occur during high heat and humidity, may bias studies based on daily data \
(@baldwin2023).

### Physiological knowledge 
Physiological evidence indicates that heat and humidity increase health risks, leading to higher rates of morbidity, mortality, reduced physical work capacity, and other adverse health outcomes as the climate crisis worsens. (@buzan2020)

The human body has two primary cooling mechanism to control the body temperature. One mechanism redistributes the blood flow to the skin (vasodilation), and therefore allows metabolically generated heat from inside of the body to reach the skin’s surface and subsequently the environment. This process requires an increased blood flow and elevates cardiac demand. This may trigger cardiovascular events such as myocardial infarction in subjects with heart conditions.
The second mechanism involves secreting sweat onto the skin. Subsequently, the sweat evaporates to the environment and heat from the body is removed. These two mentioned cooling mechanisms are crucial for the human body to maintain a sustainable core temperature. (@ebi2021)

The effectiveness of these cooling mechanisms decline with the increase of heat and humidity. A higher ambient humidity reduces the proportion of sweat that evaporates from the skin into the air. Hence, the effectiveness of sweating decreases and sweating does not lead to heat loss anymore, causing a rise in body temperature. The human body responds to such ineffective sweating by an even higher rate of sweating. This can lead to heat stroke, dehydration, or the body's physical limits may simply be exceeded. (@baldwin2023)

Certain populations are at higher risk due to less efficient cooling mechanisms. For instance, people older than 65 years may have a reduced sweating ability same as individuals with cardiovascular diseases who cannot effectively increase blood flow to the skin. These factors, among others, exacerbate the difficulty in cooling the body. (@ebi2021) \
When the body can no longer cool itself, heat stress occurs. However, these physiological insights are not consistently reflected in epidemiological studies. Several methodological and physiological factors might explain this discrepancy.

Firstly, heat-vulnerable individuals often have impaired sweating capacities, especially with aging this capacity decreases. Thermoregulatory sweating can decrease by up to 25% in individuals older than 60. Consequently, they will be less sensitive to high humidity. Since most studies only include individuals over 65 years, this might explain the lack of association between humidity and heat stress.

Secondly, it is possible that the threshold at which differences in humidity impact the body’s regulatory mechanisms are rarely attained. This means, that only at a specific threshold of humidity the body’s cooling mechanisms have trouble and this threshold may be reached rarely in especially dry areas.

Thirdly, heat-related outcomes may result from other causes unrelated to humidity. For example, cardiovascular diseases are a leading cause of death, and high temperatures can exacerbate these conditions independently of humidity.

Lastly, dehydration, a major cause of heat-related deaths, reduces the amount of sweat, thereby minimizing the impact of high humidity on evaporation and cooling. (@baldwin2023)

## Heat and Humidity in Studies 
### Humidity definitions 
Humidity has multiple definitions. This makes it crucial to select an appropriate measure to answer specific research questions. (@baldwin2023) categorize humidity variables into the four following groups: 

1. Simple relative variables including Relative Humidity (RH), Dew-point depression (DPD) and Saturation deficit. \
2. Mass-based variables including specific humidity, absolute humidity, dew point temperature, vapor pressure and mixing ratio. \
3. Composite indicators such as the Heat index, Humidex, Wet-bulb-temperature, UTCI and further more. \
4. Physiologic-based humidity indicators including maximum evaporation, required evaporation, total evaporation. 

The simple relative variable RH (%) is the ratio of actual vapor pressure (e) in the air to the saturation vapor pressure ($e_s$), indicating proximity to saturation ($RH=e/e_s$). The dew-point depression (DPD) (°C) indicates the relative dryness of the air and is the difference between the near-surface air temperature ($T_a$) and the dew point temperature ($T_d$) ($DPD=T_a-T_d$). Saturation deficit refers to the amount by which the water vapor must be increased in a given environment to reach saturation, assuming unchanged temperature and pressure.

The mass-based variable specific humidity (q) (g/kg) is the ratio of the mass of water vapor (g) to the total mass of moist air (kg). Absolute humidity is the density of water vapor, defined as the ratio of the mass of water vapor (g) to the volume of the air mixture ($m^3$). The dew point temperature ($T_d$) (°C) is the temperature at which saturation occurs. The mixing ratio (w) refers to the ratio of the mass of water vapor (g) to the mass of dry air (kg).

The Heat Index (HI) describes a "feels like" temperature, while the Humidex is a discomfort index. The wet-bulb temperature is the temperature at which air becomes saturated through evaporation at constant pressure. The Universal Thermal Climate Index (UTCI) (°C) is a physiologically based thermal stress scale designed to reflect the average human physiological response to the outdoor thermal environment.

The physiological-based humidity indicator maximum evaporation ($E_{max}$) represents the maximal evaporation rate in a given environment and clothing conditions. The required evaporation ($E_{req}$) measures the evaporation needed to achieve thermal balance, and the total evaporation (E) is the total amount of heat lost via evaporation of sweat from the body.  

Davis et al. reviewed the use of humidity variables in epidemiological studies. They found that 62% of the studies used RH, while absolute and specific humidity were used in only 5.3% and 1.5% of the studies, respectively. They discussed the appropriate use of different humidity variables, especially the use of RH due to its strong inverse correlation with temperature and its diurnal and seasonal variability was mentioned. It is recommended to use RH sparingly. Instead, they recommend using mass-based variables such as absolute or specific humidity, as they provide a more accurate measure of exposure to humidity levels. \
The selection of humidity variables must be done carefully, they can be highly correlated with other variables, for example with heat. Davis et al. concluded that using RH as humidity variable is often inappropriate and should be done with caution if necessary. (@davis2016)

The widespread use of RH in epidemiological studies is largely caused by the available data and the measurements which are performed. Mostly, if humidity is measured, then RH is measured. Deriving mass water-vapor based variables such as absolute and specific humidity can introduce uncertainty and bias in estimates.  

### Composite Indicators 
Composite indicators are frequently used in studies. There are over 120 heat stress metrics available (@buzan2020), most of which were originally designed for heat warning systems. Composite indicators consider heat, humidity, and additional variables simultaneously, making them practical for certain applications. Often they are just a weighted sum of different variables. Using such composite indices is very practical, but it comes at the cost of interpretability. Changes in the composite index can result from variations in heat, humidity, or other included variables. Therefore, information about individual effects is lost and the effects can not be separated. Furthermore, the use of such indices may lead to biased or misleading estimates. Hence, the effects of heat or humidity may be underestimated. Baldwin et al. recommend the use of alternative methods such as confounding or effect modification, to separately assess the effects of heat and humidity (@baldwin2023).  

### Confounder
A confounder is a variable that affects both the exposure of interest and the outcome. Figure 2 illustrates how humidity acts as a confounder in the relationship between heat and health outcomes. 



\begin{center}\includegraphics[width=500]{Confounder} \end{center}


**Figure 2:** Humidity as confounding effect.

Humidity has been used in many epidemiological studies as a confounder when the effect of heat on health outcomes was analyzed. Researchers aim to remove the confounding influence of humidity on the relationship between heat and health outcomes by adjusting for humidity. This means the indirect effect of humidity through heat will be eliminated. Not adjusting for confounder variables, the causal effect of heat may be overestimated. This adjustment leads to a different interpretation of the estimates, as shown in the following linear model:

$$
Ε[Y│t,h]= β_0+ β_1*t+ β_2*h
$$


Y represents the health outcome of interest, t is heat, and h humidity. In this model, $β_1$ represents the conditional total effect of heat on the health outcome at any given level of humidity. $β_2$ represents the direct effect of humidity on the outcome. That is the causal effect of humidity on the outcome when heat is held fixed and thus blocking the humidity effect on heat. 

Interpretation of estimates can however get very complex and sometimes even misleading.Consequently, the causal relation between the variables needs to be considered carefully (@baldwin2023, @westreich2013) 

One significant problem of adjusting for confounding is to potentially run into multicollinearity issues. 
If heat (the exposure of interest) and humidity (the confounder) are highly correlated and both variables are included in the same regression model multicollinearity is present. Consequently, the estimates capture similar information and therefore making it difficult to interpret individual variable effects. Additionally, the precision of the estimated coefficients may be reduced. Including further variables such as calendar variables to adjust for long-term or seasonal trends can exacerbate multicollinearity if the correlation between heat and humidity varies significantly over time. Furthermore, including other variables which are highly correlated with heat and humidity will lead to further multicollinearity issues. 

Baldwin et al. suggest using alternative approaches to significance testing that do not rely solely on p-values. Furthermore, it is crucial to carefully consider the relationship between heat and humidity and their individual effects on health outcomes to avoid misleading results. (@baldwin2023)

### Effect Modification and Interaction 
Including humidity as an effect modifier or interaction with heat is another method to test an effect of humidity in the framework of negative health outcomes and was used in some epidemiological studies (@baldwin2023). Unlike confounding, which distorts relationships and should therefore be eliminated, effect modification and interaction are desirable biases. They provide valuable insights into how variables interact and influence outcomes, revealing inherent aspects of causal relationships that need clarification. 

Often, the two terms *Effect Modification* and *Interaction* are used interchangeably. However, from a causal view these two concepts differ in a counterfactual perspective and slightly differ in their definition. The counterfactual perspective involves hypothetical scenarios, asking “What if?” questions such as, “What if heat remained the same but humidity was lower? Would the health outcome change?”. (@bours2021) Counterfactual outcomes are typically written with a superscript such as $Y^t$ and answers the question: "What would the outcome Y be if the exposure was t?". Or as another example $Y^{th}$ answers the question: "What would the outcome Y be if both exposures t and h occurred together?". As in real-world data typically not all combinations of exposures are observed in each individual we are talking about hypothetical scenarios.  
Both concepts are scale dependent and their presence depends on the scale used. Here, we use the risk difference scale for formal definitions. 

*Interaction* refers to the combined causal effect of two exposures, such as heat and humidity. *Effect modification* pertains to how the causal effect of one exposure changes across different levels of a second exposure. For example, “How does the effect of heat change across different levels of humidity?”. 
Let T refer to heat, H to humidity and Y represent the outcome of interest such as heat-related death. $Y^t$ is the counterfactual outcome under exposure t and $Y^{th}$ is the counterfactual outcome under t and h. 

Humidity (H) is **effect modifier** on the risk difference scale for the effect of heat (T) on the outcome (Y) if H is not affected by T and there are two levels of T ($t_0$,$t_1$) and two levels of H ($h_0$,$h_1$), such that:

$$
Ε[Y^{t_1}│H=h_1] - E[Y^{t_0}│H=h_1] ≠ E[Y^{t_1}|H=h_0]- E[Y^{t_0}|H=h_0]
$$

Effect modification indicates that the effect of heat on the outcome differs across different levels of humidity. By recognizing effect modification, public health interventions can better address how heat exposure interacts with humidity. This enables more targeted strategies to mitigate and adapt health risks associated with specific combinations of heat and humidity.

There is an **interaction** on the risk difference scale between heat (T) and humidity (H) on the outcome (Y) if there are two levels of T ($t_0$,$t_1$) and two levels of H ($h_0$,$h_1$), such that:
$$
Ε[Y^{t_1, h_1}] - E[Y^{t_0, h_1}] ≠ E[Y^{t_1, h_0}]- E[Y^{t_0, h_0}]
$$

Interaction requires that the joint effect of heat and humidity is different from the sum of their individual effects. If interaction is present, studying each exposure separately may not fully capture the complexity of their joint influence on the outcome. (@vanderweele2009) 

According to Baldwin and colleagues, analytical approaches to estimate interaction of effect modification are similar. They suggest that by including interaction terms or effect modifiers the main question of interest needs to be considered same as the scale of interest and the potential policy implications. (@baldwin2023)

Unfortunately, the interpretation of interaction terms and effect modifiers often leads to misinterpretation and often the scale (multiplicative or additive) also is not considered correctly. Therefore, it is crucial to understand the concepts of interaction and effect modification and how to interpret them correctly. 

### Data Limitations 
The availability and quality of weather data vary across different locations all over the world. Low-income countries, like for example most african regions, often have limited access to weather data compared to high-income countries due to limited resources for data collection. Figure 3 illustrates global sub-daily weather data from the HadISD data source. This highlights the inconsistency in data coverage. Many regions in Africa, South America, and Australia have few measurements compared to Europe or the US. Such a lack of data, especially in areas prone to extreme heat events leaves gaps in our understanding and estimation of weather-related impacts on health outcomes.  

Tropical and subtropical climates are characterized by high levels of humidity and heat. Unfortunately, these climates are particularly affected by these data limitations. Missing data in tropical and subtropical regions hinder our ability to comprehend the combined effects of heat and humidity on health outcomes. Therefore, epidemiological studies often do not include locations where humidity plays a significant role due to these data gaps.

Another challenge is the temporal resolution of weather data. Many studies rely on daily mean measurements, which may overlook fluctuations in heat and humidity throughout the day. It is known that heat and humidity fluctuate throughout the day. Using daily mean measurements may underestimate the effects of heat and humidity. Furthermore, correlations between heat and humidity can vary depending on the time of day, potentially impacting the findings of epidemiological studies. 

Furthermore, typically weather stations are situated outside urban areas, often at airports. This may not capture urban heat effects. Additionally, as most measurements are taken outdoors, indoor and outdoor variations in heat and humidity are not accounted for.

Another limitation is that typically only relative humidity (RH) is measured at weather stations. This may not provide the most accurate representation of humidity. Deriving mass-based variables like specific or absolute humidity from RH may introduce uncertainties and bias in the estimates as menioned earlier. 

In summary, the lack of weather data can result in an underestimation of the impact of humidity on health outcomes. Researchers should consider these data limitations, especially when generalizing results to different countries and weather conditions. (@baldwin2023)


\begin{center}\includegraphics[width=500]{hadisdimage1} \end{center}

**Figure 3:** Station coverage and length of record in HadISD (@raymond2024)

## Examples 
Armstrong et al. conducted a study to assess the impact of heat and humidity on daily mortality during summer months. They used a time series regression model with a distributed lag nonlinear model (DNLM) for heat, using data from 445 cities across 24 countries. Numerous different models were fitted, including different humidity levels and terms and further one model with an interaction term between heat and humidity was modeled. Model comparison using the Akaike information criterion (AIC) revealed that the best fit was achieved when humidity was included as a linear term, i.e., as confounder. One very surprisingly result was that a 23% increase in relative humidity at the 99th percentile was associated with a 1.1% decrease in mortality. This was a non-significant finding contrary to the expectation that humidity leads to an increase in mortality. Sensitivity analyses involving specific humidity and dewpoint yielded similar results. Considering humidity as a potential effect modifier, by including an interaction term, was not favored based on AIC criteria (@armstrong2019).

Baldwin et al. criticize solely rely on AIC to determine whether an interaction term should be included in the model. Consequently, as no interaction term was used in the final model, the influence remains uncertain and raises questions about the value of its exclusion. The final model used by Armstrong et al. included relative humidity as a linear term, thereby treating humidity as a confounder. As highlighted earlier, employing humidity as a confounder can introduce several challenges. 

Armstrong et al.recognized some limitations of their study. They discussed the data gaps from tropical and less developed countries, the exclusive consideration of mortality as an endpoint, and the possibility of missing the association between humidity and mortality under specific conditions (@armstrong2019).

In another epidemiological study from Schwartz et al. the effect of heat and humidity on hospital admissions for heart disease and myocardial infarction among individuals over 64 years from 12 US cities was investigated. They fitted a Poisson model for each city, incorporating regression splines to control for season and barometric pressure, and adjusting for the day of the week. Despite their broad approach, they found no evidence of an effect of humidity on hospital admissions. Possible explanations for this results are the use of relative humidity and the resulting bias. Further, considering humidity as a confounder may also explain that no effect of humidity was found. Additionally, focusing solely on individuals older than 64, who are more vulnerable to heat and sweating impairments, may have biased the study's findings. (@schwartz2004)

Barreca conducted a study that found an effect of heat and humidity on mortality rates in the US. He used data from the Global Summary of the Day (GSOD) files, which provide detailed weather station data. The study analyzed data from 1973 to 2002, calculating specific humidity using a standard meteorological formula based on mean dew point and mean station pressure.

To estimate the effects of heat and humidity, Barreca employed an ordinary-least-squares model with the following equation: 

$$
MORT_{cmy} = \sum_b \beta^b*TEMP^b_{cmy} + \sum_{b'} \alpha^{b'} * HUMID^{b'}_{cmy} + X_{cmy} \gamma + \mu_{my} + \\
\phi_{cm} + \delta_{cm} * TIME + \pi_{cm}*TIME^2 + \epsilon_{kym}
$$

In this model:

- MORT represents the monthly mortality rate (per 100,000 inhabitants) in county c, year y and calendar month m.
- TEMP is a set of temperature variables indicating the fraction of days county c experiences mean temperatures within specific 10 °F bins (e.g., 50–60 °F).
- HUMID is a set of humidity variables indicating the fraction of days county c experiences mean humidity levels within specific two-grams-of-water-vapor bins (e.g., 2–4 g/kg).
- X is a vector of controls for precipitation and temperature-humidity interaction terms.
- $\mu$ represents unrestricted time effects.
- $\phi$ represents unrestricted county by calendar month fixed effects.
- $\delta_{cm}*TIME$ and $\pi_{cm} * TIME^2$ represent unrestricted county by calendar month quadratic time trends.
Barreca clustered the standard errors by the state of residence to account for the possibility that 
$\epsilon$ is correlated within states, both across counties and over time. Additionally, the equation was weighted by the county population in 2000.
He fitted different models: one model included only temperature, and another included only humidity. Both the temperature-mortality and humidity-mortality relationships followed a similar pattern. In a subsequent model, he included both TEMP and HUMID covariates. An interesting finding from this model was that both cold temperatures and low humidity levels remained significant determinants of mortality.
Finally, he fitted a model with an interaction term between heat and humidity. The estimates from this model suggest that increasing humidity levels are more dangerous at high temperatures. 
Therefore, controlling for humidity is according to Barreca particularily important in the context of predicting distributional effects of climate change. (@barreca2012)

In a separate review, Budd examined the history and limitations of the wet-bulb globe temperature (WBGT) index, widely used since its development for the United States Army and Marine Corps in 1950. The index incorporates air temperature, mean radiant temperature, absolute humidity, and air movement as basic elements of the thermal environment. Budd identifies several limitations of the index, including its underestimation of stress in environments with restricted evaporation, challenges in interpreting its results, and limitations in the accuracy of sweat measurement due to calibration methods and instrumentation. (@budd2008)

## Discussion 
The discrepancies observed in epidemiological studies regarding the impact of humidity on health outcomes may stem from various factors. It is important to understand how to appropriately use humidity in study designs and recognize the different humidity definitions. Prioritizing mass water-vapor based variables such as specific or absolute humidity over relative humidity is a crucial step. 
As discussed earlier, humidity can be considered as a confounder, effect modifier, in an interaction term with heat or component of a composite index. Each of these considerations of humidity comes with its own set of limitations.

It is essential to provide a clear research question or hypothesis when evaluating the role of humidity in studies. Researchers must be mindful of these considerations to ensure the robustness and validity of their study results.

Efforts to improve data collection are required, especially to include regions with tropical climates and low-income countries that are often underrepresented in studies. Additionally, relying solely on vulnerable populations may introduce bias, underscoring the importance of striving for a representative sample.

Humidity represents an intriguing factor worth studying further, along with heat, especially in studies investigating health effects. Planning and executing appropriate studies to examine the effects of heat and humidity requires a nuanced understanding of humidity. It is appropriate incorporation into models, and knowledge aware of potential limitations and errors.

<!--chapter:end:10-heatandhum.Rmd-->

# Risk Projections   {#he4}

*Author: Author*

*Supervisor: Helmut Kuechenhoff*

*Degree: Bachelor*   

In many studies, the effects of climate change in the future is discussed. 
In a tutorial paper (@vicedo) , statistical methods for modelling future risks due to climate devlopement are presented. In a two other papers paper the future effects of ozone and heat  related are presented (@domingo,@chen).   



<!--chapter:end:11-projections.Rmd-->

# Open issue: Future health risks of climate change {#he5}

*Author: Author*

*Supervisor: Helmut Kuechenhoff*

*Degree: Bachelor Master  *   

Further literature should be checked. 
Possible starting point : 
Report of the IPCC- panel (@IPCC 3000 pages) 



<!--chapter:end:12-future.Rmd-->

# Open issue

*Author: Author*

*Supervisor: Helmut Kuechenhoff*

*Degree: Bachelor Master  *   

Further aspects of climate change 

Economic, political aspects, communication 

@campbell @sun 

<!--chapter:end:13-open.Rmd-->

# Acknowledgements

The most important contributions are from the students themselves.
The success of such projects highly depends on the students.
And this book is a success, so thanks a lot to all the authors!
The other important role is the supervisor.
Thanks to all the supervisors who participated!
Special thanks to [SUPERVISING PROFESSOR](https://www.stablab.stat.uni-muenchen.de/personen/leitung/kuechenhoff1/index.html) who enabled us to conduct the seminar in such an experimental way, supported us and gave valuable feedback for the seminar structure.
Thanks a lot as well to the entire [Department of Statistics](https://www.statistik.uni-muenchen.de/) and the [LMU Munich](http://www.en.uni-muenchen.de/index.html) for the infrastructure.

The authors of this work take full responsibilities for its content.



<!--chapter:end:98-acknowledgments.Rmd-->






<!--chapter:end:99-references.Rmd-->
