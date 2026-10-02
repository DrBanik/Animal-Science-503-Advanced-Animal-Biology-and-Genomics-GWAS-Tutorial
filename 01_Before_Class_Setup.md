# Before Class Setup — GWAS Workshop

Please complete the items below **before class** so we can spend the workshop on the analysis rather than installation and connection problems.

## ✅ Before class to-do list

Please try to complete these **before class** so we do not spend the whole workshop troubleshooting installations 😅

1. Download the **WSU VPN (GlobalProtect)**
2. Download **FileZilla**
3. **Windows users:** install **WSL / Ubuntu**
4. **Mac users:** install **XQuartz**
5. Contact **Mosope** if you have any questions about installation **before the class**

> ⚠️ **Important:** You must be connected to the **WSU VPN** to access the workshop server.

Your assigned username will be `student01`, `student02`, etc. Replace `studentXX` in all commands below with your assigned student number.

For example:

```bash
ssh -Y student01@10.104.58.24
```

---

## 🔐 Connect to the WSU VPN

If you are **not connected to the WSU network**, connect to the WSU VPN before attempting to log in to the workshop server.

WSU uses **GlobalProtect VPN**.

### 1. Download GlobalProtect

Open the WSU GlobalProtect download page in your web browser:

👉 https://vpn.wsu.edu/global-protect/getsoftwarepage.esp

Sign in using your **WSU Network ID** if prompted and complete multifactor authentication.

Download the appropriate **GlobalProtect** installer for your operating system.

### 🪟 Windows

Most modern Windows computers use the **64-bit** installer.

After downloading:

1. Open the downloaded installer.
2. Click **Next**.
3. Keep the default installation location unless you have a reason to change it.
4. Continue through the installation.
5. Click **Close** when installation is complete.
6. Open **GlobalProtect** from the Windows Start menu if it does not start automatically.

### 🍎 macOS

Download the **Mac GlobalProtect agent**.

After downloading:

1. Open the downloaded installer.
2. Run the **GlobalProtect Installer**.
3. Click **Continue** through the installation steps.
4. Install the GlobalProtect package.
5. If macOS asks you to approve the Palo Alto Networks system extension, open **System Settings** and allow it.
6. Complete the installation.

### 2. Connect to the WSU VPN

Open **GlobalProtect**.

When asked for the portal address, enter:

```text
vpn.wsu.edu
```

> ⚠️ Enter only `vpn.wsu.edu` in GlobalProtect. Do **not** enter `https://vpn.wsu.edu`.

Click **Connect**.

Sign in using your **WSU Network ID credentials** and complete multifactor authentication if prompted.

When the connection is successful, GlobalProtect should show:

```text
Connected
```

> 💡 You only need the VPN when you are not already connected to the WSU network.

---

## 💻 Computer setup and server login

### 🪟 Windows Users

#### 1. Update PowerShell

Please make sure that you have an updated version of PowerShell installed.

You can check your version by opening PowerShell and running:

```powershell
$PSVersionTable.PSVersion
```

If needed, download/update PowerShell here:

👉 https://learn.microsoft.com/en-us/powershell/scripting/install/installing-powershell-on-windows

#### 2. Install or update WSL

Open **PowerShell as Administrator**.

If WSL is not already installed, run:

```powershell
wsl --install
```

This installs WSL and Ubuntu by default.

If WSL is already installed, update it:

```powershell
wsl --update
```

You can check the installation with:

```powershell
wsl --version
wsl --list --verbose
```

Restart your computer if Windows asks you to do so.

#### 3. Start Ubuntu and create your Linux user

Open **Ubuntu** from the Start menu.

The first time Ubuntu starts, it will ask you to create a Linux username and password.

> 💡 This username is for your **own computer** and does not need to match your workshop `studentXX` username. Make sure you save the username and password in a secure folder.

Once Ubuntu opens, update the installed packages:

```bash
sudo apt update
sudo apt upgrade -y
```

Install the utilities needed for the workshop and X11 authentication:

```bash
sudo apt install -y xauth x11-apps openssh-client
```

#### 4. Update WSLg

From **Windows PowerShell**, not from inside Ubuntu, run:

```powershell
wsl --update
wsl --shutdown
```

Open Ubuntu again.

WSLg provides graphical Linux application support, so Windows users using current WSL generally do **not** need to install Xming or VcXsrv separately.

You can test graphical support from Ubuntu with:

```bash
xclock
```

A small clock window should appear. 🕐

#### 5. Connect to the workshop server

From the Ubuntu/WSL terminal:

```bash
ssh -Y studentXX@10.104.58.24
```

For example:

```bash
ssh -Y student01@10.104.58.24
```

Enter your assigned workshop password when prompted.

After connecting, you should see a prompt similar to:

```text
student01@gwaa-workshop:~$
```

You can confirm that you are in the correct environment with:

```bash
hostname
whoami
pwd
```

For `student01`, for example, you should see:

```text
gwaa-workshop
student01
/home/student01
```

To test X11 forwarding after connecting:

```bash
xclock
```

A clock window should appear on your Windows desktop. ✅

---

### 🖥️ PuTTY — optional for Windows users

Windows users who prefer a graphical SSH client may also install **PuTTY**.

Download PuTTY here:

👉 https://www.putty.org/

Configure the connection as:

```text
Host Name: 10.104.58.24
Port:      22
Protocol:  SSH
```

When prompted for a username, enter your assigned `studentXX` username.

> ℹ️ **Note:** PuTTY is convenient for normal terminal access, but PuTTY alone does not provide an X server. If you need graphical/X11 applications, using **Ubuntu through WSL/WSLg** is encouraged.

---

### 🍎 macOS Users

macOS already includes an SSH client, so WSL and PuTTY are **not required**.

#### 1. Install XQuartz

macOS requires an X11 server for displaying graphical applications forwarded from the workshop server.

Download and install **XQuartz** here:

👉 https://www.xquartz.org/

After installing XQuartz, **log out of macOS and log back in** or restart the computer.

Open XQuartz before testing X11 forwarding.

If you already use Homebrew, you can also install XQuartz from the command line:

```bash
brew install --cask xquartz
```

#### 2. Connect using Terminal

Open the macOS **Terminal** application and connect with:

```bash
ssh -Y studentXX@10.104.58.24
```

For example:

```bash
ssh -Y student01@10.104.58.24
```

Enter your assigned workshop password.

You should arrive at:

```text
student01@gwaa-workshop:~$
```

Test the connection:

```bash
hostname
whoami
pwd
```

Then test graphical/X11 forwarding:

```bash
xclock
```

The clock should appear through XQuartz. 🕐✅

---

### 📂 File transfer — Windows and macOS

The simplest way to transfer files is to use `scp` directly from the terminal on your local computer.

To download a **single file** from the workshop server:

```bash
scp studentXX@10.104.58.24:/path/to/file .
```

To download a whole directory, add `-r`:

```bash
scp -r studentXX@10.104.58.24:/path/to/directory .
```

The `.` at the end means download the file or directory to your current local directory.

For example:

```bash
scp student01@10.104.58.24:/workshop/students/student01/results.txt .
```

or for a whole directory:

```bash
scp -r student01@10.104.58.24:/workshop/students/student01/results .
```

💡 You can check your current local directory before downloading with:

```bash
pwd
```

Otherwise, Windows and macOS users are encouraged to download **FileZilla Client** for an easier graphical way to transfer files between your computer and the workshop server.

Download FileZilla here:

👉 https://filezilla-project.org/

In FileZilla, use:

```text
Protocol: SFTP - SSH File Transfer Protocol
Host:     10.104.58.24
Port:     22
Username: studentXX
Password: Your assigned workshop password
```

Your personal workshop directory is:

```text
/workshop/students/studentXX
```

For example:

```text
/workshop/students/student01
```

Your Linux home directory is:

```text
/home/studentXX
```

> 💡 In FileZilla, make sure the **local folder on your computer is writable** (for example, your Downloads folder) before downloading files.

---

### 🚀 Quick reference

**Windows — recommended login with X11 support:**

```bash
# Run from Ubuntu/WSL
ssh -Y studentXX@10.104.58.24
```

**macOS — recommended login with XQuartz installed:**

```bash
ssh -Y studentXX@10.104.58.24
```

**Command-line login without graphical forwarding:**

```bash
ssh studentXX@10.104.58.24
```

**FileZilla/SFTP:**

```text
Host:     10.104.58.24
Port:     22
Protocol: SFTP
Username: studentXX
```

**Your personal workshop directory:**

```text
/workshop/students/studentXX
```

✅ Please verify that you can **connect to the WSU VPN if needed, log in to the server, transfer a file using FileZilla, and open `xclock` using X11 forwarding before the workshop.**

---

