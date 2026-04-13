import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import numpy as np
from sklearn.inspection import permutation_importance
from itertools import cycle
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.model_selection import train_test_split, GridSearchCV

from imblearn.over_sampling import SMOTE
import os
traits="age"
traits2="region"

os.chdir("/Users/freya/Documents/RProjects/MSI/20250311-WT-3-age/random_forest/"+traits)
# file_path = './converted_data.xlsx'
file_path = '/Users/freya/Documents/RProjects/MSI/20250311-WT-3-age/heatmap/wide_data.xlsx'

data = pd.read_excel(file_path)
data[traits] = data[traits].astype('category')
data = data.drop(columns=[traits2])
data = data.drop(columns=["bio"])
data = data.dropna()
#删除 traits2 指定的列（如 age）


# 数据预处理
X = data.drop(traits, axis=1)
y = data[traits]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=30, stratify=y)#m6da第一，0。72

#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.6, random_state=300, stratify=y)#m6da第一，0。81
#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=150, stratify=y)#有m6da

#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=150, stratify=y)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# PCA降维
pca = PCA(n_components=3)  # 选择适当的成分数量
X_train_pca = pca.fit_transform(X_train_scaled)
X_test_pca = pca.transform(X_test_scaled)

# 训练随机森林模型
model = RandomForestClassifier(random_state=123)
model.fit(X_train_pca, y_train)
y_score = model.predict_proba(X_test_pca)

# 检查模型性能
y_pred = model.predict(X_test_pca)
print(confusion_matrix(y_test, y_pred))
print(classification_report(y_test, y_pred))
# 计算每个类别的ROC曲线并插值，使得所有的fpr和tpr具有相同的长度
from sklearn.metrics import roc_curve

fpr_all = []
tpr_all = []
roc_auc_all = []

# 选择一个固定的fpr作为基准
mean_fpr = np.linspace(0, 1, 100)

for i, class_name in enumerate(y.cat.categories):
    fpr, tpr, _ = roc_curve(y_test == class_name, y_score[:, i])
    roc_auc = auc(fpr, tpr)

    # 插值，确保所有的fpr和tpr长度一致
    tpr_interp = np.interp(mean_fpr, fpr, tpr)
    fpr_all.append(mean_fpr)
    tpr_all.append(tpr_interp)
    roc_auc_all.append(roc_auc)

# 计算平均tpr
tpr_avg = np.mean(tpr_all, axis=0)
roc_auc_avg = auc(mean_fpr, tpr_avg)

# 绘制合并后的ROC曲线
plt.figure(figsize=(4, 4))
plt.plot(mean_fpr, tpr_avg, color='blue', lw=2, label=f' (AUC = {roc_auc_avg:.2f})')

# 绘制虚线（随机猜测线）
plt.plot([0, 1], [0, 1], 'k--', lw=2)

# 设置坐标轴和标签
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontsize=14)
plt.ylabel('True Positive Rate', fontsize=14)
plt.title('Mean ROC Curve ('+traits+')', fontsize=14)
plt.legend(loc='lower right', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)

# 保存ROC曲线图
roc_curve_path = 'mean_roc_curve.pdf'
plt.tight_layout()
plt.savefig(roc_curve_path, format='pdf', bbox_inches='tight')
plt.show()

# 打印保存的ROC曲线图路径
print(f'Mean ROC curve saved to {roc_curve_path}')
# 计算所有类别的特征重要性并合并为一个图
all_importances = []

# 计算每个类别的permutation importance并合并
pca_components = pd.DataFrame(pca.components_, columns=X.columns, index=[f'PC{i + 1}' for i in range(pca.n_components_)])

for class_name in y_test.cat.categories:
    result = permutation_importance(model, X_test_pca, (y_test == class_name).astype(int), n_repeats=10, random_state=123, n_jobs=2)
    class_importances = result.importances_mean

    if np.all(class_importances == 0):
        print(f"No significant features found for class: {class_name}")
        continue  # 跳过没有显著特征的类别

    # 计算原始特征的重要性
    original_feature_importance = np.abs(pca_components.T.dot(class_importances))

    # 将每个类别的特征重要性添加到all_importances列表中
    all_importances.append(original_feature_importance)



# 计算特征重要性
feature_importances = model.feature_importances_
# 获取PCA成分在原始特征中的权重
pca_components = pd.DataFrame(pca.components_, columns=X.columns, index=[f'PC{i+1}' for i in range(pca.n_components_)])
# 计算每个原始特征的综合重要性
original_feature_importance = np.abs(pca_components.T.dot(feature_importances))
# 创建一个DataFrame来存储原始特征及其重要性
importance_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': original_feature_importance
})
# 排序
importance_df = importance_df.sort_values(by='Importance', ascending=False)
# 选择前10个最重要的特征
top_10_importance_df = importance_df.head(10)
# 绘制重要性图
plt.figure(figsize=(5, 4))  # 调整为更小的图形

plt.barh(top_10_importance_df['Feature'], top_10_importance_df['Importance'], color='skyblue', height=0.6)
plt.xlabel('Importance', fontsize=12)  # 增大x轴标签字体
plt.ylabel('Features', fontsize=12)   # 增大y轴标签字体
plt.title('Feature importances ('+traits+')', fontsize=14)  # 增大标题字体
plt.gca().invert_yaxis()

# 设置行标签的字体大小
plt.yticks(fontsize=10)  # 增大y轴的字体大小
plt.xticks(fontsize=10)  # 增大x轴的字体大小

# 保存图像
feature_importance_path = 'feature_importance.pdf'
plt.savefig(feature_importance_path, format='pdf', bbox_inches='tight')
plt.show()
plt.close()

