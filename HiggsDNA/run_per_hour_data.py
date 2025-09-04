#!/usr/bin/env python3
import subprocess
import time
import os
import signal
from datetime import datetime

# ===== 參數 =====
CMD = ["bash", "scripts/run_analysis_data_run2.sh"]
# NUM_CYCLES = 5            # 週期次數
# NUM_CYCLES = 10           # 週期次數
NUM_CYCLES = 10            # 週期次數
# WORK_SECONDS = 3600       # 單輪跑多久（秒）= 1 小時
WORK_SECONDS = 1800       # 單輪跑多久（秒）= 30 分鐘
# WORK_SECONDS = 900        # 單輪跑多久（秒）= 15 分鐘
# WORK_SECONDS = 600        # 單輪跑多久（秒）= 10 分鐘
RESTART_DELAY = 10        # 殺掉後等多久再重啟（秒）
TERM_GRACE = 20           # 給 SIGTERM 的寬限（秒），之後升級 SIGKILL

def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def start_proc():
    # 在新 session 啟動（等同新 process group），方便一次殺整組
    proc = subprocess.Popen(
        CMD,
        preexec_fn=os.setsid,  # Unix：建立新 session / process group
        stdout=None, stderr=None, stdin=None,
    )
    print(f"[{now()}] 啟動：PID={proc.pid}  指令={' '.join(CMD)}")
    return proc

def kill_proc_group(proc):
    if proc.poll() is not None:
        print(f"[{now()}] 行程已結束（returncode={proc.returncode}），略過終止。")
        return

    pgid = os.getpgid(proc.pid)
    try:
        print(f"[{now()}] 發送 SIGTERM 給 Process Group PGID={pgid}…")
        os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError:
        print(f"[{now()}] 找不到行程群組（可能已退出），略過。")
        return

    # 給一段寬限時間
    waited = 0
    while waited < TERM_GRACE:
        if proc.poll() is not None:
            print(f"[{now()}] 子行程在 {waited} 秒內正常結束（returncode={proc.returncode}）。")
            return
        time.sleep(1)
        waited += 1

    if proc.poll() is None:
        print(f"[{now()}] 仍未結束，發送 SIGKILL 給 PGID={pgid}。")
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass

        # 再等一下確認退出
        for _ in range(5):
            if proc.poll() is not None:
                break
            time.sleep(1)

    print(f"[{now()}] 最終狀態：returncode={proc.returncode}")

def main():
    print(f"[{now()}] === 開始循環（總共 {NUM_CYCLES} 次）===")
    proc = None
    try:
        for i in range(1, NUM_CYCLES + 1):
            print(f"[{now()}] --- 週期 {i}/{NUM_CYCLES} ---")
            proc = start_proc()

            # 讓它跑指定時間；若提早結束就不再等待
            t0 = time.time()
            while time.time() - t0 < WORK_SECONDS:
                if proc.poll() is not None:
                    print(f"[{now()}] 子行程提前結束（returncode={proc.returncode}），結束等待。")
                    break
                time.sleep(2)

            # 無論是否提前結束，都進入「終止 → 延遲 → 重啟」流程
            kill_proc_group(proc)
            print(f"[{now()}] 等待 {RESTART_DELAY} 秒後重啟…")
            time.sleep(RESTART_DELAY)

        print(f"[{now()}] === 所有週期完成 ===")

    except KeyboardInterrupt:
        print(f"\n[{now()}] 收到 KeyboardInterrupt，開始清理…")
        if proc is not None:
            kill_proc_group(proc)
        print(f"[{now()}] 已清理並退出。")

if __name__ == "__main__":
    main()
