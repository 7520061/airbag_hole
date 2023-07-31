function para = Para_sort(para)

vent_num = length(para);

% for j = 1 : 1 : vent_num/2
%     para(j) = para(j) * 1000;
%     para(j) = round(para(j),TieBreaker="even");
%     if rem(para(j), 2) == 1
%         para(j) = para(j) + 1;
%     end
%     para(j) = para(j) / 1000;
%     if para(j) < 0.006
%         para(j) = 0.000;
%     end
% end

for j = vent_num/2 : 1 : vent_num
    para(j) = para(j) * 1000;
    para(j) = round(para(j),TieBreaker="even");
    if rem(para(j), 2) == 1
        para(j) = para(j) + 1;
    end
    para(j) = para(j) / 1000;
end

%% 距離を降順に設計変数をソート

vent_D = (para(1:vent_num/2))';
vent_H = (para(vent_num/2+1:end))';

[vent_H_sorted,I] = sort(vent_H,'descend');
vent_D_sorted = vent_D(I);

para = [vent_D_sorted;vent_H_sorted];
para = para';

end