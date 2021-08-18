- In destpts_1101_pld60, I have FEWER particles than the ones I released. What is happening here?
	- it looks like an error associated just with destpts_1101_pld60
	- some of the points don't even look like they were tracked. Luckily most of these were from MPAs that I will exclude, but there are some significant MPAs where this is not the case (1 big one in Alaska and a few on the central coast).
	- However, the connectivity lines seem correct. There are connections FROM those MPAs that don't see to have points. This makes me think that the dataframe was correct and something happened with outputting the points.
	- Further investigating: in the results, there is no shapefile for section 5 of 1101 pld 60


SO, this is a biology script to get JUST that shapefile. It is configured for 1 nc file and 1 pld.