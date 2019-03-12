function DesignController
给定期望剪切频率
给定期望闭环幅值
给定期望闭环相位
给定双十检查频点
while 1
    n = 1
    给定期望损失相角
    给定期望剪切频率
    while问题无解
        补相位
        设置非线性环节
        if 求解成功 
            break
        else
            增加期望损失相角
            增加期望剪切频率
            if 增加到顶
                break
            end
        end
    end
    if 求解成功
        break
    else
        放宽期望闭环幅值
        放宽期望闭环相位
    end
end

            