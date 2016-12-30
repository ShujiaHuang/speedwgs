# 使用  RAM 子账号授权驻云上传 ODPS 数据的方法

##   开通 RAM 服务

访问 https://buy.aliyun.com/ram 开通 RAM 服务

##  创建 RAM 账号

在 RAM 控制台 https://ram.console.aliyun.com/ 用户管理页面点击新建用户

注意：
 - 勾选『为该用户自动生成AccessKey』
 - 点开 AccessKey 详情并及时记录生成的 Access ID/Key（如图）
![Screenshot](file:///home/lyman/create_ram.png)

## 将 RAM 账号加入 ODPS Project 并授权

RAM 账号在 ODPS 中的表达为『RAM$主账号:RAM账号名』。

本文档中新建的 RAM 账号名叫 `huada_zhuyun`，由 `odpsdemo_admin@aliyun-test.com` 主账号创建，那么在  ODPS 中描述这个子账号时应为 `RAM$odpsdemo_admin@aliyun-test.com:huada_zhuyun`

用主账号的身份运行 odpscmd，并执行以下命令来完成授权：

```
use huada_test;
add accountprovider RAM;
add user RAM$odpsdemo_admin@aliyun-test.com:huada_zhuyun;
grant CreateInstance, CreateTable, List on project huada_test to user RAM$odpsdemo_admin@aliyun-test.com:huada_zhuyun;
```

注意：
 - 替换 add user 和 grant 语句的 RAM 账号（`RAM$odpsdemo_admin@aliyun-test.com:huada_zhuyun`）
 - 替换 use 和 grant 语句的 Project（huada_test）

## 将 RAM 账号的 access id/key 交给驻云的同学

按类似以下格式编辑配置文件，填入 RAM 账号的 Access ID 和 Access Key 及 Project 名，并发送给驻云的同学，完成上传

```
access_id=xxxxxxxxxxxx
access_key=xxxxxxxxxxxxxxxxxxxxxxxxx
project_name=huada_test
end_point=http://service.odps.aliyun.com/api
tunnel_endpoint=http://dt.nu16.odps.aliyun.com
https_check=true
log_view_host=http://logview.odps.aliyun.com
```

## 销毁 RAM 账号

上传完毕后，在 RAM 控制台删除用过的账号 https://ram.console.aliyun.com
