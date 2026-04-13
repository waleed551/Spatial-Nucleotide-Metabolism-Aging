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

os.chdir("/Users/freya/Documents/RProjects/MSI/CDF/test/randomforest")
file_path = '/Users/freya/Documents/RProjects/MSI/CDF/test/randomforest/WT-10-combined.xlsx'
data = pd.read_excel(file_path)
data['sheet'] = data['sheet'].astype('category')
data = data.dropna()
#数据预处理
X = data.drop('sheet', axis=1)
y = data['sheet']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=123, stratify=y)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 使用SMOTE处理数据不平衡问题
#smote = SMOTE(random_state=123)
#X_train_res, y_train_res = smote.fit_resample(X_train_scaled, y_train)

#PCA降维
pca = PCA(n_components=3)  # 选择适当的成分数量
X_train_pca = pca.fit_transform(X_train_scaled)
X_test_pca = pca.transform(X_test_scaled)
#训练随机森林模型
model = RandomForestClassifier(random_state=123)
model.fit(X_train_pca, y_train)
y_score = model.predict_proba(X_test_pca)

# 使用网格搜索优化随机森林模型参数
#param_grid = {
#    'n_estimators': [100, 200, 300],
#    'max_depth': [None, 10, 20, 30],
 #   'min_samples_split': [2, 5, 10]
#}
#grid_search = GridSearchCV(RandomForestClassifier(random_state=123), param_grid, cv=5, n_jobs=-1, scoring='accuracy')
#grid_search.fit(X_train_pca, y_train_res)

# 最优模型
#best_model = grid_search.best_estimator_

# 检查模型性能
y_pred = model.predict(X_test_pca)
print(confusion_matrix(y_test, y_pred))
print(classification_report(y_test, y_pred))

#计算ROC曲线和AUC
# 计算每个类别的ROC曲线
fpr = dict()
tpr = dict()
roc_auc = dict()
n_classes = len(y.cat.categories)
for i, class_name in enumerate(y.cat.categories):
    fpr[class_name], tpr[class_name], _ = roc_curve(y_test == class_name, y_score[:, i])
    roc_auc[class_name] = auc(fpr[class_name], tpr[class_name])

# 绘制每个类别的ROC曲线
plt.figure(figsize=(12, 8))
colors = cycle(['aqua', 'darkorange', 'cornflowerblue', 'red', 'green', 'blue', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'magenta', 'yellow', 'black'])

for (class_name, color) in zip(y.cat.categories, colors):
    plt.plot(fpr[class_name], tpr[class_name], color=color, lw=2, label=f'{class_name} (AUC = {roc_auc[class_name]:.2f})')

plt.plot([0, 1], [0, 1], 'k--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontsize=14)
plt.ylabel('True Positive Rate', fontsize=14)
plt.title('ROC Curve for Random Forest Model(Positive and Negative)', fontsize=16)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)

# 调整布局以确保图例不重叠
plt.tight_layout(rect=[0, 0, 0.8, 1])

# 保存图像
roc_curve_path = 'roc_curve.pdf'
plt.savefig(roc_curve_path, format='pdf', bbox_inches='tight')
plt.show()

# 打印保存的ROC曲线图路径
print(f'ROC curve saved to {roc_curve_path}')


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

# 绘制重要性图
plt.figure(figsize=(14, 10))
plt.barh(importance_df['Feature'], importance_df['Importance'], color='skyblue')
plt.xlabel('Importance')
plt.ylabel('Features')
plt.title('Feature Importances(Positive and Negative)')
plt.gca().invert_yaxis()

# 设置行标签的字体大小
plt.yticks(fontsize=8)

# 保存图像
feature_importance_path = 'feature_importance.pdf'
plt.savefig(feature_importance_path, format='pdf', bbox_inches='tight')
plt.close()



#计算每个类别的贡献度
# 创建保存图片的文件夹
output_folder = './feature_importance_plots'
os.makedirs(output_folder, exist_ok=True)
# 对每个类别分别计算特征重要性
pca_components = pd.DataFrame(pca.components_, columns=X.columns, index=[f'PC{i + 1}' for i in range(pca.n_components_)])

for class_name in y_test.cat.categories:
    class_index = list(y_test.cat.categories).index(class_name)

    # 为每个类别计算 permutation_importance
    result = permutation_importance(model, X_test_pca, (y_test == class_name).astype(int), n_repeats=10, random_state=123, n_jobs=2)
    class_importances = result.importances_mean

    if np.all(class_importances == 0):
        print(f"No significant features found for class: {class_name}")
        continue

    # 计算原始特征的重要性
    original_feature_importance = np.abs(pca_components.T.dot(class_importances))

    importance_df = pd.DataFrame({
        'Feature': X.columns,
        'Importance': original_feature_importance
    })

    importance_df = importance_df.sort_values(by='Importance', ascending=False)

    plt.figure(figsize=(12, 8))
    plt.barh(importance_df['Feature'], importance_df['Importance'], color='skyblue', edgecolor='black')
    plt.xlabel('Importance', fontsize=14)
    plt.ylabel('Features', fontsize=14)
    plt.title(f'Feature Importances for class: {class_name}', fontsize=16)
    plt.gca().invert_yaxis()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=8)  # 调整行标签字体大小

    # 调整边距
    plt.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.1)

    # 保存为PDF文件
    output_path = os.path.join(output_folder, f'feature_importance_{class_name}.pdf')
    plt.savefig(output_path, format='pdf', bbox_inches='tight')
    plt.close()

# 打印生成的文件夹路径
print(f'Feature importance plots saved to {output_folder}')