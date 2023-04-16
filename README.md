## WASH Cholera

This repo contains scripts to complete the analyses of the [Water, Sanitation and Cholera in sub-Saharan Africa](https://duckduckgo.com) manuscript

**Description**

Univariate country-level and district-level analyses: `country_district_univariate.R`

Multivariate regression models: `incidence_rf_gbm.R`

Multivariate classification models: `hotspot_rf_gbm.R` 

Multivariate classification models with hotspot data: `hotspot_rf_gbm_high_incidence.R` 
 

**Country Level Data Summary**

Aggregated population weighted water and sanitation access/reliance percentage and incidence of cholera per 1,000 by country using the estimates from [ref](https://doi.org/10.1016/S2214-109X(20)30278-3) and [ref](https://doi.org/10.1016/S0140-6736(17)33050-7) for 2010-2016.  
 
<small>
| Country |Piped Water |Other Improved Water |Piped Sanitation | Other Improved Sanitation|Unimproved water|Surface Water|Unimproved Sanitation|Open Defecation|Incidence/1,000|
|----|-----|-----|----|-----|----|-----|-----|----|------|
Angola|32|26.6|32.6|9.5|11.5|29.9|28.4|29.5|0.09904
Burundi|36.4|46.4|5|47.6|11.4|5.8|43.1|4.3|0.2468
Benin|39.1|33.2|3.6|23.9|21.6|6.1|7.7|64.8|0.04665
Burkina Faso|18|58.4|1.3|29.6|21.7|1.9|3.3|65.8|0.00358
Central African Rep.|12.1|46.4|0.4|32.3|37.6|4|39.1|28.2|0.02886
CÃ´te d'Ivoire|42.6|35.2|15.8|29.2|14.3|7.8|17.2|37.7|0.01891
Cameroon|35.7|34.6|9.7|44.6|20.7|9|37.2|8.5|0.50913
DR Congo|21.2|22.3|3.9|35.9|45.3|11.2|43.9|16.3|0.36399
Republic of Congo|27.1|34.2|4.8|26.2|21.3|17.4|50.3|18.8|0.20317
Ethiopia|37.5|23.2|1.9|17.2|23.4|16|40.6|40.3|0.06755
Gabon|61.9|20.7|20.9|16|1.8|15.6|49.3|13.8|0.00033
Ghana|28.3|54.9|13.3|48.3|4.7|12.1|7.6|30.9|0.40156
Guinea|23.9|49.8|9.5|32.3|13|13.3|38.6|19.6|0.20086
Gambia|63|23.4|7.3|56.2|13.4|0.2|12.3|24.1|0.00612
Guinea-Bissau|23.8|34.7|7|46|41.1|0.4|20.8|26.1|0.39827
Equatorial Guinea|39|30.8|8.3|33.8|16.3|13.8|50.9|7|0.001
Kenya|31.5|28.5|6.2|41.5|13.4|26.6|33.4|18.9|0.17645
Liberia|2.9|69.8|9.8|22.8|7.5|19.8|12.5|54.9|0.11201
Madagascar|23.9|19.1|3.2|8.6|32.7|24.3|36.8|51.4|0.00895
Mali|32.6|34.4|4.1|37.8|29.1|3.9|44.1|14|0.0326
Mozambique|39.9|12.5|1.3|24.2|33|14.6|35.1|39.4|0.12538
Mauritania|48.1|31|7.1|38.4|20.4|0.5|11.6|42.9|0.01136
Malawi|21.3|62.1|1.5|40.9|12.6|4.1|47.7|9.9|0.03755
Namibia|73.2|8.5|35.9|7.6|10.5|7.7|3.4|53.1|0.17558
Niger|39.6|28|4.7|18.4|31.4|1|14.1|62.8|0.06494
Nigeria|11.1|53.1|15.8|27.1|19.4|16.4|22|35.2|0.10247
Rwanda|33.9|39.8|0.6|84.1|15.1|11.2|13.4|1.9|0.01835
Sudan|28.3|44.7|0.1|20.3|21.5|5.6|27|52.7|0.00598
Senegal|56.8|12.3|15.8|53.7|30.5|0.4|6.9|23.6|0.00072
Sierra Leone|16.1|40.3|4.8|42.7|19.4|24.2|28|24.5|1.89451
Somalia|32.3|52.8|7.6|27.5|5.3|9.6|7.8|57.1|1.27452
South Sudan|3.2|30.9|0.7|12.1|48.7|17.2|18.9|68.4|0.61066
Swaziland|54.3|16.9|16.2|62.3|10.1|18.7|8.7|12.8|0.00065
Chad|18.9|36.3|1.6|13.5|37.5|7.2|15.4|69.5|0.86626
Togo|27.6|33.9|10.4|20.3|18.9|19.7|9.2|60.1|0.01796
Tanzania|38.2|20.4|4.9|29.1|25|16.4|49|17|0.18485
Uganda|18|57.8|0.2|78.7|13.6|10.6|10.2|10.9|0.03673
South Africa|87.6|3.9|55.4|19.8|2.3|6.2|20.8|4|0.00103
Zambia|28.4|28.9|3.7|32.2|27.4|15.4|43.6|20.6|0.11792
Zimbabwe|32.9|42.6|32.3|29.4|14.6|9.9|6.9|31.4|0.03265

