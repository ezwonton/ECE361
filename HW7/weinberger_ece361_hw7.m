
receiver1=[4.324792635
1.886096562
1.318833992
2.610091288
2.42254717
1.804664798
1.158303063
0.878588478
0.914926586
2.39231983
1.557308839
1.82536144
3.28445573
2.298036125
1.630987157
1.198121862
3.235369279
2.517372324
3.741562988
2.936084419
1.589694429
1.390686808
1.10270262
1.110041901
2.631286934
2.183205305
2.326764688
3.183181795
4.073587036
1.990822039
2.749405032
0.712195023
3.009151424
3.019081308
3.772506241
2.570356796
2.13221093
2.022436004
1.731047395
1.483172022
1.410103227
2.41231523
1.826828271
1.368774606
1.684376342
1.455473755
1.627343368
1.388252745
0.885738513
0.801420428
1.058130083
1.656284269
0.593032474
2.712827219
1.201362759
2.068218952
1.584072632
1.08810766
1.840045611
3.112815605
1.214151795
0.685936004
1.250657973
2.61797645
2.377232109
1.645084817
2.664314676
2.310647943
1.481726305
0.713202181
2.593836338
1.638987436
1.401492496
1.915216565
2.185810742
1.077236448
1.851871399
1.718702696
2.142120147
1.426259113
2.157675358
2.262592712
2.465338792
2.558557544
2.814900609
3.660748324
3.096748915
1.849987894
2.518208296
2.999521187
0.702715164
2.865756046
2.063289286
1.992487107
1.215953197
1.408585535
2.864983849
0.902170848
2.841277925
1.785965633
1.639614176
1.944195363
1.938886246
1.529896028
2.394710374
1.929458314
0.643344392
2.551186877
1.477136425
0.225580038
1.295162753
1.766588325
3.506613733
2.30167346
2.269637225
3.150698798
0.750179077
0.639571481
2.039236456
2.005706438
1.005044307
0.553111812
2.076822277
1.515623613
2.064233361
1.496322112



0.959492026
3.562337081
2.389676287
0.814954166
3.227697776
1.534479176
0.554061663
3.401310469
1.061249835
1.128169052
1.781685356
0.476561091
1.616983644
2.315627641
1.910938685
1.256497773
0.608356495
2.652981766
1.723692855
1.702357819
1.640342379
1.974166279
1.049247975
1.701092085
0.60065944
1.239264617
2.008611626
0.917308626
0.626606131
0.479094165
1.390235104
2.446342408
1.556175284
1.833746086
3.819410297
1.929943862
1.527685335
1.295587383
1.82613363
1.639374446
2.445877615
3.110274509
1.944180016
4.013385786
1.8618902
1.488786528
2.596267387
0.448438248
2.133359247
2.015361794
2.5539819
1.273445148
1.015239021
1.851734635
2.922210884
2.014213873
1.612103677
1.467826952
2.941159665
1.449124121
0.660276846
1.962699917
2.825607161
1.457340125
1.58625822
1.50525168
1.59641103
1.842632661
1.760012532
1.598384178
1.817568618
3.531800595
0.568411632
1.756171148];
[data_ks, data_points]=ksdensity(receiver1,'NumPoints', 200);
plot(data_points, data_ks)
ndata=mean(receiver1)/std(receiver1);

hold on
arithmean=zeros(199,1);
for i=1:199
    arithmean(i,1)=(receiver1(i,1)+receiver1(i+1))/2;
end
[arithmean_ks, arithmean_points]=ksdensity(arithmean,'NumPoints', 199);
plot(arithmean_points, arithmean_ks)
narith=mean(arithmean)/std(arithmean);

hold on
geomean=zeros(199,1);
for i=1:199
    geomean(i,1)=sqrt(receiver1(i,1)*receiver1(i+1));
end
[geomean_ks, geomean_points]=ksdensity(geomean,'NumPoints', 199);
plot(geomean_points, geomean_ks)
ngeo=mean(geomean)/std(geomean);
hold on
maxx=zeros(199,1);
for i=1:199
    maxx(i,1)=max(receiver1(i,1),receiver1(i+1));
end
[maxx_ks, maxx_points]=ksdensity(maxx,'NumPoints', 199);
plot(maxx_points, maxx_ks)
nmaxx=mean(maxx)/std(maxx);
hold on
xlabel('values')
ylabel('estimated density')
title('data (weinberger)')
words1=sprintf('input:mean/std. dev=%.2f',ndata);
words2=sprintf('AM:mean/std. dev =%.2f',narith);
words3=sprintf('GM:mean/std. dev =%.2f',ngeo);
words4=sprintf('MAX:mean/std. dev =%.2f',nmaxx);
legend(words1,words2,words3,words4)